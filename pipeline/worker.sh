#!/bin/bash
./StereoGene cfg=stereogene.cfg 1.lst
./Confounder cfg=stereogene.cfg 1.lst
./Projector cfg=confounder.cfg 1.lst
./StereoGene cfg=stereogene1.cfg 1.lst






