#! /usr/bin/env Rscript --arch=x86_64

library(RUnit)
source("/Users/brad/Rcode/fmrireg/pkg/workspace.R")


testsuite <- defineTestSuite("fmrireg", dirs = "tests",
	testFileRegexp = "^runit.+\\.R",
	testFuncRegexp = "^test.+",
	rngKind = "Marsaglia-Multicarry",
	rngNormalKind = "Kinderman-Ramage")

testResult <- runTestSuite(testsuite)
printTextProtocol(testResult)
