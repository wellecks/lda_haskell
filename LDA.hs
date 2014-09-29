{- LDA.hs

Sean Welleck

Online Latent Dirichlet Allocation.
-}

{-# LANGUAGE FlexibleContexts #-}

module LDA where

import Data.Matrix
import qualified Data.Vector as DV (Vector, toList, length, fromList) 
import qualified Data.IntMap as IM
import Control.Monad.State

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
-- A mapping from wordID -> wordCount
type WordCounts = IM.IntMap Int

-- | A Document
--    * ws the raw text
--    * wCts a mapping of wordID -> wordCount for the document
data Document = Document {
  ws   :: [String],
  wCts :: WordCounts
}

data Batch = Batch {
  docs    :: [(Int, Document)]
}

-- Aliases for matrices.
type Lambda   = Matrix Double
type ELogBeta = Matrix Double
type SStats   = Matrix Double
type Gamma    = Matrix Double

-- ================= MODEL TRAINING ================= 

-- Update a given model using the given batch
update :: Model -> Batch -> RandomModel Model
update m batch = do
  m' <- eStep batch
  _  <- put m'
  let tau0   = tau $ hyperParams m
  let count0 = updateCount m
  let kappa0 = kappa $ hyperParams m
  let rho'      = rhot tau0 count0 kappa0
  let lambda'   = updateLambda m (sstats m') (d m')
  let eLogBeta' = dirichletExpectation lambda'
  return  Model {
            k = k m,
            w = w m,
            d = d m,
            updateCount = updateCount m + 1,
            gamma       = gamma m',
            eLogBeta    = eLogBeta',
            expELogBeta = mexp eLogBeta',
            lambda      = lambda',
            vocab       = vocab m,
            sstats      = sstats m',
            hyperParams = HyperParams {
              alpha = alpha $ hyperParams m,
              eta   = eta   $ hyperParams m,
              rho   = rho',
              kappa = kappa $ hyperParams m,
              tau   = tau   $ hyperParams m
              }
          }

-- Lambda update from the paper.
updateLambda :: Model -> SStats -> Int -> Lambda
updateLambda m ss batchSize = elementwise (+) scaleLambda0 scaleSS
  where lambda0 = lambda m 
        r = rho $ hyperParams m
        e = eta $ hyperParams m 
        nd = fromIntegral $ d m 
        bs = fromIntegral batchSize
        scaleLambda0 = scaleMatrix (1.0 - r) lambda0
        scaleSS      = scaleMatrix r $ madd e (scaleMatrix (nd / bs) ss)

-- TODO
-- Expectation step
eStep :: Batch -> RandomModel Model
eStep batch = do
  m <- get
  (gamma', sstats') <- eStepOuter batch
  let sstats'' = elementwise (*) sstats' (expELogBeta m)
  return  Model {
            k = k m,
            w = w m,
            d = d m,
            updateCount = updateCount m,
            gamma       = gamma',
            eLogBeta    = eLogBeta m,
            expELogBeta = expELogBeta m,
            lambda      = lambda m,
            vocab       = vocab m,
            sstats      = sstats'',
            hyperParams = HyperParams {
              alpha = alpha $ hyperParams m,
              eta   = eta $ hyperParams m,
              rho   = rho $ hyperParams m,
              kappa = kappa $ hyperParams m,
              tau   = tau $ hyperParams m
              }
          }

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
  let alph     = alpha $ hyperParams m  
  let beta     = expELogBeta m
  return $ foldl (eStepDocIter alph eeltheta beta) (gamma0, sstats0) (docs batch)

-- Perform e-step update for a single document in a batch.
eStepDocIter :: Double -> Matrix Double -> Matrix Double -> (Gamma, SStats) -> (Int, Document) -> (Gamma, SStats)
eStepDocIter alph expELogTheta expeLogBeta (gamm, ss) (i, doc) = 
  let ids = IM.keys $ wCts doc
      cts = DV.fromList $ map fromIntegral $ intMapVals $ wCts doc 
      gammad        = getRow i gamm
      expElogthetad = getRow i expELogTheta
      expElogbetad  = selectCols ids expeLogBeta
      phinorm   = updatePhinorm expElogthetad expElogbetad
      newgammad = innerIter (100::Int) gammad alph expElogthetad cts phinorm expElogbetad
      newsstats = updateSStats ss ids expElogthetad cts phinorm
  in  (updateRow i gamm newgammad, newsstats)   

innerIter :: Int -> DV.Vector Double -> Double -> DV.Vector Double -> DV.Vector Double ->  DV.Vector Double -> Matrix Double -> DV.Vector Double
innerIter n g a e c p b = 
  if n <= 0 
  then g
  else innerIter (n-1) (updateGamma g a e c p b) a e c (updatePhinorm e b) b

-- Note that since RVar is the inner monad, we lift the RVar computation into
-- the StateT monad.
randomGammaMatrixM :: Int -> Int -> RandomModel (Matrix Double)
randomGammaMatrixM r c = do
  xs <- lift $ replicateM (r*c) $ Data.Random.sample gammaGen
  return $ fromList r c xs

-- sstats[:, ids] += n.outer(expElogthetad.T, cts/phinorm)
updateSStats :: SStats -> [Int] -> DV.Vector Double -> DV.Vector Double ->  DV.Vector Double -> SStats
updateSStats ss ids elthetad c p = updateCols ids1 (multStd2 (transpose tmat) (elementwise (/) cmat pmat)) ss
  where cmat = fromList 1 (DV.length c) $ DV.toList c
        pmat = fromList 1 (DV.length p) $ DV.toList p
        tmat = fromList 1 (DV.length elthetad) $ DV.toList elthetad
        ids1 = map (+1) ids

updateGamma :: DV.Vector Double -> Double -> DV.Vector Double -> DV.Vector Double ->  DV.Vector Double -> Matrix Double -> DV.Vector Double
updateGamma _ a e c p b = getRow 1 $ madd a $ elementwise (*) emat $ multStd2 (elementwise (/) cmat pmat) (transpose b)
  where cmat = fromList 1 (DV.length c) $ DV.toList c
        pmat = fromList 1 (DV.length p) $ DV.toList p
        emat = fromList 1 (DV.length e) $ DV.toList e

updatePhinorm :: DV.Vector Double -> Matrix Double -> DV.Vector Double
updatePhinorm e m = getRow 1 $ multStd2 emat m
  where emat = fromList 1 (DV.length e) $ DV.toList e

-- ========= DATA LOADING and DATA FORMATTING =======

-- TODO
docToWordCounts :: Document -> WordCounts
docToWordCounts = undefined

-- print the topics for a model
printTopics :: Model -> [String]
printTopics m = [show (updateCount m), show (lambda m)]

-- load the vocabulary file
loadVocabulary :: String -> [String]
loadVocabulary = undefined

docListFromFile :: String -> [Document]
docListFromFile = undefined

-- TODO
makeBatches :: [Document] -> Vocabulary -> [Batch]
makeBatches ds _ = [Batch { docs = dtuples }]
  where dtuples = zip (take (length ds) (iterate (+1) 1)) ds

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
initModel :: Lambda -> Model
initModel lambda0 = 
  let k' = 10
      w' = 16
      d' = 4 in
  Model {
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
  lambda0 <- randomGammaMatrix 10 16
  let model0  = initModel lambda0
  let batches = makeBatches sampleDocs vocabulary
  finalModel <- Data.Random.runRVar (evalStateT (foldM update model0 batches) model0) Data.Random.StdRandom
  print $ printTopics finalModel
