
module MathUtils where

import Data.Matrix
import qualified Data.Vector as DV (Vector, toList, length, foldr, fromList, zip) 
import qualified Data.IntMap as IM
import Data.List (sort, (\\))
import Control.Monad.State
import Numeric.Digamma

import Data.Random

-- ================= MATRICES ==================== 

-- Convert a vector to a 1 x n matrix.
vecMat :: DV.Vector a -> Matrix a
vecMat v = fromList 1 (DV.length v) $ DV.toList v

-- Add a scalar to every element of a matrix.
madd :: Double -> Matrix Double -> Matrix Double
madd x = mapAll $ \y -> y + x

-- e^ each element of the matrix.
mexp :: Matrix Double -> Matrix Double
mexp = mapAll $ \x -> exp 1 ** x

-- Map over every element of a matrix.
mapAll :: (a -> b) -> Matrix a -> Matrix b
mapAll f m = fromList (nrows m) (ncols m) $ map f $ Data.Matrix.toList m 

intMapVals :: IM.IntMap a -> [a]
intMapVals im = map (\x -> im IM.! x) $ IM.keys im

-- Sums the rows of a Matrix.
sumRows :: Num a => Matrix a -> [a]
sumRows m = map (sum . rowList) [1 .. nrows m]
  where rowList i' = DV.toList $ getRow i' m

-- Set row i in the matrix to the values in the provided vector.
-- Note that Data.Matrix matrices are 1-indexed.
updateRow :: Int -> Matrix a -> DV.Vector a -> Matrix a
updateRow i m v = foldl (fun i) m (zip [1 .. ncols m] (DV.toList v))
  where fun :: Int -> Matrix a -> (Int, a) -> Matrix a
        fun r mat (c, val) = setElem val (r, c) mat

-- Update the specified columns using the values in the given matrix.
updateCols :: Show a => [Int] -> Matrix a -> Matrix a -> Matrix a
updateCols idx um m = foldr (updateCol um) m $ zip [1 .. (ncols um-1)] (sort idx)

-- Set column mInd of m to column umInd of um.
-- Concatenate m and um and swaps columns, then remove um.
updateCol :: Matrix a -> (Int, Int) -> Matrix a -> Matrix a
updateCol um (umInd, mInd) m = condense $ swapCols $ m <|> um
  where condense = submatrix 1 (nrows m) 1 (ncols m)
        swapCols = switchCols mInd (ncols m + umInd)

-- Mimics numpy matrix[:, idx], where idx is a list of indices.
selectCols :: (Fractional a) => [Int] -> Matrix a -> Matrix a
selectCols idx m = foldr removeCol m $ [1 .. ncols m] \\ sort idx

-- Mimics numpy matrix[:, ~idx], where idx is a list of indices.
filterCols :: (Fractional a) => [Int] -> Matrix a -> Matrix a
filterCols idx m = foldr removeCol m $ sort idx

removeCol :: (Fractional a) => Int -> Matrix a -> Matrix a
removeCol n m = m --fromList (nrows m) (ncols m - 1) $ DV.foldr (\(x, i) acc -> if i /= n && (n - (i `mod` (ncols m))/= 0) then x:acc else acc) [] $ DV.zip (getMatrixAsVector m) $ DV.fromList [1..(ncols m * nrows m)]

-- minorMatrix (nrows m + 2) n $ setSize 0.0 (nrows m + 1) (ncols m) m
  -- where sub1 = submatrix 1 (nrows m) 1 (n-1) m
  --       sub2 = submatrix 1 (nrows m) (n+1) (ncols m) m

-- Subtract a column vector from each column.
subFromCol :: Num a => [a] -> Matrix a -> Matrix a
subFromCol col m = elementwise (-) m $ repeatCol col (ncols m)

-- Repeat a column vector n times to create a matrix
repeatCol :: [a] -> Int -> Matrix a
repeatCol col n = horizN initial initial (n - 1)
  where initial = fromList (length col) 1 col 
        horizN :: Matrix a -> Matrix a -> Int -> Matrix a
        horizN x _ 0 = x
        horizN x y m = horizN (x <|> y) y (m - 1)

-- ================= MATH / PROB. ==================== 

-- How much to weight the information from a mini-batch.
-- p_t = (t_0 + t)^{-k}
rhot :: Double -> Int -> Double -> Double
rhot tau0 count kap = (tau0 + fromIntegral count) ** (-kap)

-- For a vector theta ~ Dir(alpha), computes E[log(theta)] given alpha.
dirichletExpectation :: Matrix Double -> Matrix Double
dirichletExpectation m = mapAll digamma $ subFromCol (sumRows m) m

zeroDouble :: Int -> Int -> Matrix Double
zeroDouble = zero

gammaGen :: RVar Double
gammaGen = gamma 100.0 (1.0/100.0)

randomGammaMatrix :: (MonadRandom m) => Int -> Int -> m (Matrix Double)
randomGammaMatrix r c = do
  xs  <- replicateM (r*c) $ sample gammaGen
  return $ fromList r c xs