# We declare the variables to be visible inside the function only in the 
# function signature. These variables have underscore as prefix, come after
# the user variables and they are meant to be ignored by the user.  

function getSamTag(tag, _idxkey, _key){
    if(NF < 12){
        return ""
    }
    for(i=12; i <= NF; i++){
        _idxkey= index($i, ":")
        _key=substr($i, 1, _idxkey-1)
        if(_key == tag){
            return substr($i, length(_key) + 4)
        }
    }
}

function getInfoTag(tag, value_idx, _target, _retval, _eqIdx, _key, _list){
	if(NF < 8){
		return 0
	}
    split($8, _target, ";")
    _retval= 0
    for(i in _target){
        _eqIdx= index(_target[i], "=")
        if(_eqIdx == 0){
            if(tag == _target[i]){
                _retval= 1 # VCF tag of type boolean: No associated value
                break
            }
        } 
        else {
            _key= substr(_target[i], 1, _eqIdx-1)
            if(_key == tag){
                _retval= substr(_target[i], _eqIdx+1)
                break 
            }
        }
    }
    if(value_idx == ""){
        return _retval
    } else {
    	split(_retval, _list, ",")
    	return _list[value_idx]
    }
}

function getFmtTag(tag, sample_idx, value_idx, _fmt_array, _tagIdx, _sample_data, _list){
    # tag:        Tag to extract
    # sample_idx: Return the tag value of this sample index
    if(sample_idx == ""){
    	sample_idx= 1
    }
    if(NF < sample_idx+9 || sample_idx < 1 || sample_idx !~ /^[0-9]+$/){
        return "" # Not enough VCF fields/wrong input
    }
    split($9, _fmt_array, ":")
    _tagIdx= -1
    for(i in _fmt_array){
        if(_fmt_array[i] == tag){
            _tagIdx= i
            break
        }
    }
    if(_tagIdx == -1){
        return ""
    }
    split($(9+sample_idx), _sample_data, ":")
    if(value_idx == ""){
        return _sample_data[_tagIdx]	
    } else {
    	split(_sample_data[_tagIdx], _list, ",")
    	return _list[value_idx]
    }
}