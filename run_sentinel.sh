#!/bin/sh

cd /home/user/sentinel

echo ".First <- function(){
file.remove('.Rprofile');
library(utils); library(stats); 
library(graphics); library(grDevices); 

source('include.R');
showMeWarnings(source('sentinel.R')); 

cat('\nPress ENTER to close window...\n\n'); readline();
q(save='no')}" > .Rprofile
R
rm .Rprofile
