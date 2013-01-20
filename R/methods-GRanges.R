rOverlaps <- function(gr1, gr2, fr){
	o1 <- findOverlaps(gr1, gr2)
	i <- queryHits(o1)
	j <- subjectHits(o1)
	if(length(i) > 1){
		g1 <- gr1[i, ]
		g2 <- gr2[j, ]
		ui <- unique(i)
		## some ranges not in gr1
		if(length(ui) < length(gr1)){
			gr1.only <- gr1[-i, ]
		}
		x <- intersect(g1, g2)
		n.probes.intersection <- countOverlaps(x, fr)
		p1 <- n.probes.intersection/numberProbes(g1)
		p2 <- n.probes.intersection/numberProbes(g2)
		ans <- sum(p1 > 0.5 & p2 > 0.5)/length(gr1)
	} else {
		ans <- 0
	}
	return(ans)
}
reciprocalOverlap <- function(gr1, gr2, fr){
	gr1 <- gr1[numberProbes(gr1) >= 10 & state(gr1) %in% c(1,2,5,6), ]
	gr2 <- gr2[numberProbes(gr2) >= 10 & state(gr2) %in% c(1,2,5,6), ]
	p1 <- rOverlaps(gr1, gr2, fr)
	p2 <- rOverlaps(gr2, gr1, fr)
	c(p1, p2)
}

usePennCNVTrioState <- function(x){
	f <- as.integer(substr(x, 1,1))
	m <- as.integer(substr(x, 2,2))
	o <- as.integer(substr(x, 3,3))
	indexScale <- function(x){
		y <- x
		if(any(y <= 2)) {
			x[y<=2] <- y[y<=2]+1
		}
		if(any(y > 2)){
			x[y>2] <- y[y>2]+2
		}
		return(x)
	}
	f <- indexScale(f)
	m <- indexScale(m)
	o <- indexScale(o)
	paste(f,m,o, sep="")
}
