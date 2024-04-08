#!/usr/bin/env bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

# download the code package
# commitHash="4bc2badcd31401bcdea8be5a2efc778bdf99fc57" # for xscape 1.1

# using a commit from the iSS repository that is compatible with X-SCAPE 1.1.1
folderName="iSS"
commitHash="37ce2ae5c6c5a1d7c4d762e42fabab5bc450cbe8"

rm -fr iSS
git init $folderName
cd $folderName
git remote add origin https://github.com/USPHydro/iSS.git
git fetch --depth=1 origin $commitHash
git checkout $commitHash
rm -fr $folderName/.git
