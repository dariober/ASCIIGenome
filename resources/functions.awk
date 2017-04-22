function getSamTag(tag){
    if(NF < 12){
        return ""
    }
    for(i=12; i <= NF; i++){
        idxkey= index($i, ":")
        key=substr($i, 1, idxkey-1)
        if(key == tag){
            return substr($i, length(key) + 4)
        }
    }
}