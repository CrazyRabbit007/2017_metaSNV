{
    if( length(_[FNR]) < 1) {
        _[FNR]=_[FNR] sprintf("%s\t",$1)
    }    

    _[FNR]=_[FNR] sprintf("%s\t",$2)
}

END{
    for(x=1; x < ARGC; x++) {
	split(ARGV[x],a,".cov")
	printf("\t%s",a[1])
    }
    printf("\n")
    for(x=1; x < NR/(ARGC-1)+1; x++) {
        print _[x]
    }
}