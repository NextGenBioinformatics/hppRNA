


library(Rsubread)

rm(list=ls(all=TRUE))

args <- commandArgs(TRUE);

buildindex(basename=args[1],reference=args[2]);


