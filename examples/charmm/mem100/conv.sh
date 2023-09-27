#!/bin/bash

  grep CRY test.pdb > mem.pdb
  grep POPC test.pdb >> mem.pdb
  echo "TER" >> mem.pdb
  grep DOPC test.pdb >> mem.pdb
  echo "TER" >> mem.pdb
  grep TIP3 test.pdb >> mem.pdb
  echo "TER" >> mem.pdb
  echo "END" >> mem.pdb
