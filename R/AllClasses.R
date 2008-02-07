setClassUnion("NullOrNumeric", c("numeric", "NULL"))
setClassUnion("NullOrNumericOrMatrix", c("NULL", "numeric", "matrix"))
setClassUnion("NullOrMatrix", c("NULL", "matrix"))
setClassUnion("CharacterOrNumeric", c("character", "numeric"))

setClass("HmmOptions", representation(SnpClass="character",
                                      states="character",
                                      cn.location="numeric",
                                      cn.robustSE="CharacterOrNumeric",
                                      cn.ICE="logical",
                                      gt.ICE="logical",
                                      gt.gte="CharacterOrNumeric",
                                      gte.state="CharacterOrNumeric",
                                      gt.confidence.states="character",
                                      annotation="character",
                                      notes="character"))

setClass("HmmParameter", representation(hmmOptions="HmmOptions",
                                        tau="matrix",  ##transition probability
                                        tau.scale="matrix",
                                        pi="numeric",  ##initialStateProbability
                                        beta="matrix", ##emission probability
                                        featureData="AnnotatedDataFrame",
                                        notes="character"))

setClass("HmmPredict", representation(states="character",
                                      predictions="matrix",
                                      breakpoints="list",
                                      SnpClass="character",
                                      featureData="AnnotatedDataFrame"))

##Classes for graphical parameters
setClass("ParHSet", contains="ParESet")
setClass("ParHmmSnpCallSet", contains="ParHSet")
setClass("ParHmmSnpCopyNumberSet", contains="ParHSet")
setClass("ParHmmSnpSet", contains="ParHSet")

###########################################################################
##Deprecated classes
###########################################################################
##setClass("hSet", representation(states="character",
##                                predictions="matrix",
##                                breakpoints="list"))
##setClass("HmmSnpSet", contains="hSet")







                        
         
           
           
           
                        
                        
