setMethod("calculateEmission", "oligoSnpSet", function(object,
						       mu,
						       probs,
						       states,
						       probMissing,
						       envir,
						       verbose=FALSE,
						       ice=FALSE){
	if(all(is.na(cnConfidence(object)))){
		sds <- robustSds(copyNumber(object))
	} else {
		sds <- 1/cnConfidence(object)
	}
	emission.cn <- copynumberEmission(copynumber=copyNumber(object),
					  mu=mu,
					  sds=sds,
					  states=states,
					  takeLog=FALSE,
					  verbose=verbose)
	if(!ice){
		emission.gt <- genotypeEmission(genotypes=calls(object),
						probHomCall=probs,
						probMissing=probMissing,
						states=states)
	} else {
		##assumed order
		## ROH, normal
		emission.gt <- genotypeEmissionCrlmm(genotypes=calls(object),
						     conf=callsConfidence(object),
						     pHetCalledHom=probs[1],
						     pHetCalledHet=probs[2],
						     pHomInNormal=probs[3],
						     pHomInRoh=probs[4],
						     annotation=annotation(object))
		emit.gt <- array(NA, dim=dim(emission.cn))
		dimnames(emit.gt) <- dimnames(emission.cn)
		emit.gt[, , c(1, 3)] <- emission.gt[, , "ROH"]
		emit.gt[, , c(2, 4)] <- emission.gt[, , "norm"]
		emission.gt <- emit.gt
	}
	envir[["emission.cn"]] <- emission.cn
	envir[["emission.gt"]] <- emission.gt
})

setMethod("calculateEmission", "SnpCopyNumberSet", function(object, mu, states, envir, ...){
	if(all(is.na(cnConfidence(object)))){
		sds <- robustSds(copyNumber(object))
	} else {
		sds <- 1/cnConfidence(object)
	}
	emission.cn <- copynumberEmission(copynumber=copyNumber(object),
					  mu=mu,
					  sds=sds,
					  states=states,
					  takeLog=FALSE,
					  verbose=verbose)
	envir[["emission.cn"]] <- emission.cn
})

setMethod("calculateEmission", "SnpCallSet", function(object, probs, probMissing, states, envir, ice=FALSE){
	if(!ice){
		emission.gt <- genotypeEmission(genotypes=calls(object),
						probHomCall=probs,
						probMissing=probMissing,
						states=states)
	} else{
		emission.gt <- genotypeEmissionCrlmm(genotypes=calls(object),
						     conf=callsConfidence(object),
						     pHetCalledHom=probs["pHetCalledHom"],
						     pHetCalledHet=probs["pHetCalledHet"],
						     pHomInNormal=probs["pHomInNormal"],
						     pHomInRoh=probs["pHomInRoh"])
	}
	envir[["emission.gt"]] <- emission.gt
})






				    

