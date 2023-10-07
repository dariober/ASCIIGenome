# We declare the variables to be visible inside the function only in the 
# function signature. These variables have underscore as prefix, come after
# the user variables and they are meant to be ignored by the user.  

function _trim(_x, _char){
	# For iternal use only.
    # Return a string where leading and trailing char have been removed from
    # string x. If _char is missing trim whitespaces.
    if(_char == ""){
        _char= "[ \t]"
    }
    gsub("^"_char"+", "", _x)
    gsub(_char"+$", "", _x)
    return _x
}

## Functions for SAM records
## =========================

function getSamTag(tag,    _idxkey, _key){
	# If there are no tags the sam record has 11 fields. However, we edited the raw sam record 
	# so instead of testing for 12 fields we account for the extra ones. 
	# See Track.prepareSAMRecordForAwk
    # Tags start at 12th column
    for(i=12; i <= NF; i++){
        _idxkey= index($i, ":")
        _key=substr($i, 1, _idxkey-1)
        if(_key == tag){
            return substr($i, length(_key) + 4)
        }
    }
    return ""
}

function bitset(intquery,    intqry, samflag, qrybit, sambit){
    intqry= intquery
    samflag= $2
    if(samflag !~ "^[0-9]+$" || intqry !~ "^[0-9]+$"){
        return 0    
    }
    while(intqry > 0){
        qrybit= intqry % 2
        intqry -= qrybit
        intqry = intqry / 2
        sambit= samflag % 2
        if(samflag > 0){
            samflag -= sambit
            samflag = samflag / 2
        }
        if(qrybit == 1 && sambit == 0){
            return 0
        }
    }
    return 1
}

function nobitset(intquery,    intqry, samflag, qrybit, sambit){
    intqry= intquery
    samflag= $2
    if(samflag !~ "^[0-9]+$" || intqry !~ "^[0-9]+$"){
        return 0    
    }
}

function getAlnEnd(_alnEnd) {
	_alnEnd = getSamTag("$alnEnd") # Tag must must match Track.prepareSAMRecordForAwk
	return _alnEnd
}

# Pure awk implementation is ~ 3x slower on Jawk
#function getAlnEnd(_alnEnd, _value, _cigar) {
#    if($6 == "*") {
#        return -1
#    }
#    _alnEnd = $4 - 1
#    _value = ""
#    split($6, _cigar, "")
#    for (i=1; i <= length($6); i++) {
#        if(_cigar[i] ~ /^[0-9]$/){
#            _value = _value _cigar[i]
#        } else {
#            if("MDN=X" ~ _cigar[i]) {
#                _alnEnd += _value
#            }
#            _value=""
#        }
#    }
#    return _alnEnd
#}

function getAlnLen(_alnLen) {
    if($6 == "*") {
        return -1
    }
    _alnLen = (getAlnEnd() - $4) + 1
    return _alnLen
}


## Functions for VCF
## =================

function getInfoTag(tag, value_idx,    _target, _retval, _eqIdx, _key, _list){
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

function getFmtTag(tag, sample_idx, value_idx,    _fmt_array, _tagIdx, _sample_data, _list){
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

## Functions for GFF/GTF/BED

function getGtfTag(tag,    _attrs, _attr, _i, _tagval, _n, _retval){

    if(NF < 9){
        return "" # Not enough fields
    }
    _attrs= _trim($9, ";")
    split(_attrs, _attr, ";")

	for(_i in _attr){
        _attr[_i]=  _trim(_attr[_i])
        _n= split(_attr[_i], _tagval, " ")
            
        _tagval[1]= _trim(_tagval[1])
        if(_tagval[1] == tag){
            _retval= _trim(_trim(_tagval[_n]), "\"")
            break
        }
    }
    return _retval
}

function getGffTag(tag, value_idx,    _attrs, _attr, _tagval, _i, _n, _vals, _retval){
    if(NF < 9){
        return "" # Not enough fields
    }
    _attrs= _trim($9, ";")
    split(_attrs, _attr, ";")
    for(_i in _attr){
        _attr[_i]=  _trim(_attr[_i])
        _n= split(_attr[_i], _tagval, "=")

        _tagval[1]= _trim(_tagval[1])
        if(_tagval[1] == tag){
            if(value_idx > 0){
                # Return only one of the comma-separated values
                split(_trim(_tagval[_n]), _vals, ",")
                _retval= _trim(_vals[value_idx])
            } else {
                _retval= _trim(_tagval[_n])
            }
            break
        }
    }
    return _retval
}