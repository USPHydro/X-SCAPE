#!/usr/bin/env bash
commitHash=4becb350456ae33721c8c001ca7203969ecb3cfb
folderName=AMPTGenesis
rm -fr $folderName
git clone git@github.com:USPHydro/AMPTGenesis.git --depth=1 -b $commitHash $folderName
