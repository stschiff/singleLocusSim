{-# LANGUAGE OverloadedStrings #-}

import Control.Monad.Trans.State.Strict (State, get, put, evalState)
import Control.Monad (replicateM)
import System.Random (StdGen, random, mkStdGen)
import Turtle


data CmdOpts = CmdOpts {
    optInitialAlleleCount :: Int,
    optPopSize :: Int,
    optSelection :: Double,
    optIterations :: Int,
    optSeed :: Int
} deriving (Show)

cmdOptParser :: Parser CmdOpts
cmdOptParser = CmdOpts <$> parseInitialAlleleCount <*> parsePopSize <*>
        parseSelection <*> parseIterations <*> parseSeed
  where
    parseInitialAlleleCount = optInt "initialCount" 'x' "the initial allele count of the mutation (default=1)" <|> pure 1
    parsePopSize = optInt "popSize" 'N' "the effective haploid population size (default=200)" <|> pure 200
    parseSelection = optDouble "selection" 's' "the selection coefficient on the mutation"
    parseIterations = optInt "iterations" 'i' "the number of iterations"
    parseSeed = optInt "seed" 'S' "random seed"


main :: IO ()
main = do
    CmdOpts x0 n s nIt seed <- options "single locus simulation tool"
        cmdOptParser
    let results = evalState (replicateM nIt (runSingleIt n s x0))
            (mkStdGen seed)
    let fixations = map fst results
        nFixations = sum $ map fromEnum fixations
        fixTimes = [t | (fix, t) <- results, fix]
        extinctTimes = [t | (fix, t) <- results, not fix]
        fixationProb = fromIntegral nFixations / fromIntegral nIt
        meanFixTime = fromIntegral (sum fixTimes) /
            fromIntegral (max nFixations 1)
        meanExtTime = fromIntegral (sum extinctTimes) /
            fromIntegral (max (nIt - nFixations) 1)
    echo . unsafeTextToLine $ format ("fixationProb\t"%g) fixationProb
    echo . unsafeTextToLine $ format ("meanFixTime\t"%g) meanFixTime
    echo . unsafeTextToLine $ format ("meanExtTime\t"%g) meanExtTime

runSingleIt :: Int -> Double -> Int -> State StdGen (Bool, Int)
runSingleIt n s x0 = go 0 x0
  where
    go gen x | x == n = return (True, gen)
             | x == 0 = return (False, gen)
             | otherwise = do
                 x' <- nextGen n s x
                 go (gen + 1) x'

nextGen :: Int -> Double -> Int -> State StdGen Int
nextGen n s x = do
    let freq = fromIntegral x / fromIntegral n
    getBinomialRandom (freq * exp s) n

getBinomialRandom :: Double -> Int -> State StdGen Int
getBinomialRandom prob n = do
    bernoulliResults <- replicateM n (getBernoulliRandom prob)
    return . sum . map fromEnum $ bernoulliResults

getBernoulliRandom :: Double -> State StdGen Bool
getBernoulliRandom prob = do
    stdGen <- get
    let (rn, stdGen') = random stdGen
    put stdGen'
    return $ rn < prob
