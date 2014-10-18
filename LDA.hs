{- LDA.hs

Sean Welleck

Online Latent Dirichlet Allocation.
-}

module LDA where

import Data.Matrix
import qualified Data.Vector as DV (Vector, fromList, map) 
import qualified Data.IntMap as IM
import qualified Data.Map as Map
import Data.List.Split (chunksOf)
import Data.Char (isAlphaNum)
import Control.Monad.State.Strict
import System.IO
import qualified Data.Random

import MathUtils

-- ==================== TYPES =======================

-- | The state of the LDA model at a point in time.
--    * hyperParams
--    * k           number of topics
--    * w           number of words in the vocabulary
--    * d           number of documents
--    * updateCount cumulative number of updates
--    * gamma       parameterizes theta; size d x k 
--    * eLogBeta    E[log(Beta)], where Beta is a matrix of p(w|topic)
--    * expELogBeta e^{E[log(Beta)]}
--    * lambda      parameterizes Beta
--    * vocab       the vocabulary
data Model = Model {
  hyperParams  :: HyperParams,
  k            :: Int,
  w            :: Int,
  d            :: Int,
  updateCount  :: Int,
  gamma        :: Matrix Double,
  eLogBeta     :: Matrix Double, 
  expELogBeta  :: Matrix Double,
  lambda       :: Matrix Double,
  sstats       :: SStats,
  vocab        :: Vocabulary
} deriving (Show)

-- | Hyperparameters for an LDA model.
--    * alpha Prior on topic weights Theta
--    * eta   E[log(Beta)], where Beta is a matrix of p(w|topic)
--    * rho 
--    * kappa learning parameter; decay factor for influence of batches
--    * tau   learning parameter to downweight early documents
data HyperParams = HyperParams {
  alpha :: Double,
  eta   :: Double,
  rho   :: Double,
  kappa :: Double,
  tau   :: Double
} deriving (Show)

-- A monadic model that can keep state and have random variables.
type RandomModel a = StateT Model Data.Random.RVar a
-- A mapping from wordID -> word
type Vocabulary = IM.IntMap String
-- A mapping from word   -> wordID
type Dictionary = Map.Map String Int
-- A mapping from wordID -> wordCount
type WordCounts = IM.IntMap Int

-- | A Document
--    * ws the raw text
--    * wCts a mapping of wordID -> wordCount for the document
data Document = Document {
  ws   :: [String],
  wCts :: WordCounts
} deriving (Show)

data Batch = Batch {
  docs :: [(Int, Document)]
}

type Filename = String

-- Aliases for matrices.
type Lambda   = Matrix Double
type ELogBeta = Matrix Double
type SStats   = Matrix Double
type Gamma    = Matrix Double

-- ================= MODEL TRAINING ================= 

-- Update a given model using the given batch.
update :: Model -> Batch -> RandomModel Model
update _ batch = do
  eStep batch
  m' <- get
  let m2 = updateLambda m' $ length $ docs batch
  let m3 = updateBeta m2
  let m4 = updateRho m3
  let m5 = incrModel m4
  return m5

-- Lambda update from the paper.
updateLambda :: Model -> Int -> Model
updateLambda m batchSize = do
  --m <- get
  let lambda0 = lambda m 
  let r = rho $ hyperParams m
  let e = eta $ hyperParams m 
  let nd = fromIntegral $ d m 
  let bs = fromIntegral batchSize
  let scaleLambda0 = scaleMatrix (1.0 - r) lambda0
  let scaleSS      = scaleMatrix r $ madd e (scaleMatrix (nd / bs) (sstats m))
  let lambda'      = elementwise (+) scaleLambda0 scaleSS
  m { lambda = lambda' }

updateBeta :: Model -> Model --RandomModel ()
updateBeta m = do
  -- m <- get
  let elb' = dirichletExpectation $ lambda m
  m { eLogBeta = elb', expELogBeta = mexp elb' }

updateRho :: Model -> Model -- RandomModel ()
updateRho m = do
  --m <- get
  let rho' = rhot (tau $ hyperParams m) (updateCount m) (kappa $ hyperParams m)
  m { hyperParams = (hyperParams m) { rho = rho' }}

incrModel :: Model -> Model -- RandomModel ()
incrModel m = do
  --m <- get
  m { updateCount = updateCount m + 1 }

-- Expectation step. Updates gamma and sstats.
eStep :: Batch -> RandomModel ()
eStep batch = do
  m       <- get
  (g, ss) <- eStepOuter batch
  let ss' =  elementwise (*) ss (expELogBeta m)
  put $ m { gamma = g, sstats = ss' }

-- In the e-step double loop, we're ultimately computing values for
--      gamma
--      sstats
eStepOuter :: Batch -> RandomModel (Gamma, SStats)
eStepOuter batch = do
  m <- get
  gamma0 <- randomGammaMatrixM (d m) (k m)
  let eltheta  = dirichletExpectation gamma0
  let eeltheta = mexp eltheta
  let sstats0  = zeroDouble (nrows (lambda m)) (ncols (lambda m))
  foldM (eStepDocIter eeltheta) (gamma0, sstats0) (docs batch)

-- Perform e-step update for a single document in a batch.
-- * et expELogTheta
eStepDocIter :: Matrix Double -> (Gamma, SStats) -> (Int, Document) -> 
                RandomModel (Gamma, SStats)
eStepDocIter et (g, ss) (i, doc) = do
  m <- get
  let alph = alpha $ hyperParams m
  let ids = IM.keys $ wCts doc
  let cts = DV.fromList $ map fromIntegral $ intMapVals $ wCts doc 
  let gammad  = getRow i g
  let etd     = rowVector $ getRow i et
  let ebd     = selectCols ids $ expELogBeta m
  let phinorm = updatePhinorm etd ebd
  let gammad' = innerIter 100 gammad alph etd cts phinorm ebd
  let ss'     = updateSStats ss ids etd cts phinorm
  return (updateRow i g gammad', ss')   

-- * n; g gammad; a alpha; e expElogthetad; c cts; p phinorm; b expElogbetad
innerIter :: Int -> DV.Vector Double -> Double -> Matrix Double -> 
             DV.Vector Double -> Matrix Double -> Matrix Double -> 
             DV.Vector Double
innerIter n g a e c p b = if n <= 0 then g 
  else innerIter (n-1) (updateGamma a e c p b) a e c (updatePhinorm e b) b

-- sstats[:, ids] += n.outer(expElogthetad.T, cts/phinorm)
updateSStats :: SStats -> [Int] -> Matrix Double -> DV.Vector Double ->  
                Matrix Double -> SStats
updateSStats ss ids elthetad c p = updateCols ids1 dot ss
  where tmat  = elthetad
        ids1  = map (+1) ids
        cDivP = elementwise (/) (vecMat c) p
        dot   = multStd2 (transpose tmat) cDivP

updateGamma :: Double -> Matrix Double -> 
               DV.Vector Double ->  Matrix Double -> Matrix Double -> 
               DV.Vector Double
updateGamma a e c p b = toVec $ madd a eMultDot
  where emat = e
        dot  = multStd2 cDivP (transpose b)
        eMultDot = elementwise (*) emat dot
        cDivP    = elementwise (/) (vecMat c) p
        toVec    = getRow 1

updatePhinorm :: Matrix Double -> Matrix Double -> Matrix Double
updatePhinorm e m = madd 1e-50  $ multStrassenMixed e m

randomGammaMatrixM :: Int -> Int -> RandomModel (Matrix Double)
randomGammaMatrixM r c = do
  xs <- lift $ replicateM (r*c) $ Data.Random.sample gammaGen
  return $ fromList r c xs

-- ========= DATA LOADING and DATA FORMATTING =======

countWords :: Dictionary -> [String] -> WordCounts
countWords dict = foldr (\s m -> IM.insertWith (+) (dict Map.! s) 1 m) IM.empty

-- Print the topics for a model
printTopics :: Model -> [String]
printTopics m = [show (updateCount m), show (lambda m)]

-- Produce a vocabulary from the given file.
loadVocabulary :: Filename -> IO Vocabulary
loadVocabulary filename = do
  ls <- readWholeFile filename
  return $ IM.fromList $ zip [0..] ls

-- Produce a dictionary from the given file.
loadDictionary :: Filename -> IO Dictionary
loadDictionary filename = do
  ls <- readWholeFile filename
  return $ Map.fromList $ zip ls [0..]

readWholeFile :: Filename -> IO [String]
readWholeFile f = do
  s <- readFile f
  return $ lines s

readLazily :: Filename -> IO [String]
readLazily f = do
  inHandle <- openFile f ReadMode
  liftM lines $ hGetContents inHandle

splitLine :: String -> [String]
splitLine = map (filter isAlphaNum) . words

-- Read a file, creating a document from each line.
-- We need the Vocabulary to retrieve the wordIDs.
docListFromFile :: Filename -> Dictionary -> IO [Document]
docListFromFile f dict = do
  ls <- readWholeFile f
  let filtered = map (rmNonDict dict) ls
  return $ map (
    (\line -> Document { ws = line, wCts = countWords dict line}) . splitLine) filtered

rmNonDict :: Dictionary -> String -> String
rmNonDict dict l= unwords $ filter (`Map.member` dict) (splitLine l)

loadInput :: Filename -> Filename -> IO ([Document], Vocabulary)
loadInput dataFile dictFile = do
  dict <- loadDictionary dictFile
  voc  <- loadVocabulary dictFile
  ds   <- docListFromFile dataFile dict
  let nonEmpty = filter (not . null . ws) ds
  return (nonEmpty, voc)

makeBatches :: Int -> [Document] -> [Batch]
makeBatches batchSize ds = 
  map (\c -> Batch { docs = dtuples c }) $ chunksOf batchSize ds
  where dtuples = zip [1..]

sampleDocs :: [Document]
sampleDocs = [
  Document { ws = ["the", "tree", "runs", "quickly"], wCts = IM.fromList [(0,1), (1, 1), (2, 1), (3, 1)]},
  Document { ws = ["far", "away", "lions", "yell"],   wCts = IM.fromList [(4,1), (5, 1), (6, 1), (7, 1)]},
  Document { ws = ["she", "shoots", "shapes", "now"], wCts = IM.fromList [(8,1), (9, 1), (10, 1), (11, 1)]},
  Document { ws = ["he", "hits", "hops", "happily"],  wCts = IM.fromList [(12,1), (13, 1), (14, 1), (15, 1)]}]

vocabulary :: Vocabulary
vocabulary = IM.fromList [(0, "the"), (1, "tree"), (2, "runs"), (3, "quickly"), 
                          (4, "far"), (5, "away"), (6, "lions"), (7, "yell"), 
                          (8, "she"), (9, "shoots"), (10, "shapes"), (11, "now"), 
                          (12, "he"), (13, "hits"), (14, "hops"), (15, "happily")]

-- ============== MAIN EXECUTION ====================

-- create a base initial model
initModel :: Int -> Int -> Int -> Lambda -> Model
initModel k' w' d' lambda0 = Model {
    k            = k',
    w            = w',
    d            = d',
    updateCount  = 0,
    gamma        = zero d' k',
    eLogBeta     = zero k' w', 
    expELogBeta  = zero k' w',
    lambda       = lambda0,
    vocab        = vocabulary,
    sstats       = zero k' w',
    hyperParams = HyperParams {
      alpha = 1.0 / fromIntegral k',
      eta   = 1.0 / fromIntegral k',
      rho   = 0.7,
      kappa = 0.7,
      tau   = 1024
    }
  }

main :: IO ()
main = do
  let k' = 10
  (ds, voc)  <- loadInput "data/nytimes_tweets.txt" "data/dictnostops.txt"
  print $ length ds
  let b' = 50
  let d' = length ds
  let w' = IM.size voc
  lambda0    <- randomGammaMatrix k' w'
  let model0  = initModel k' w' b' lambda0
  let batches = makeBatches b' ds
  finalModel <- Data.Random.runRVar (evalStateT (foldM update model0 batches) model0) Data.Random.StdRandom
  print $ printTopics finalModel

