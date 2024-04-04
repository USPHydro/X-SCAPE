#!/usr/bin/env bash
commitHash=4becb350456ae33721c8c001ca7203969ecb3cfb
folderName=AMPTGenesis
rm -fr $folderName
git init $folderName
cd $folderName
git remote add origin git@github.com:USPHydro/AMPTGenesis.git
git fetch --depth=1 origin $commitHash
git checkout $commitHash
