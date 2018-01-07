## Simple help functions for testing bash scripts. 

assertFileExists(){ 
    # USAGE:
    # assertFileExists myfile.txt
    if [ -e "$1" ]
    then
        echo -e "\033[32mPass\033[0m"
        return 0
    else
        echo -e "\033[31mFAILED. File does not exist: $1\033[0m"
    fi
    return 1
}

assertFileIsEmpty(){
    # USAGE
    # assertFileIsEmpty err.log
    if [ ! -e "$1" ]
    then
        echo -e "\033[31mFAILED. File does not exist: $1\033[0m"
        return 1
    fi
    if [ -s "$1" ]
    then
        echo -e "\033[31mFAILED. File is not empty: $1\033[0m"
        return 1
    else
        echo -e "\033[32mPass\033[0m"
        return 0
    fi
}


assertFileIsNotEmpty(){
    # USAGE
    # assertFileIsNotEmpty err.log
    if [ ! -e "$1" ]
    then
        echo -e "\033[31mFAILED. File does not exist: $1\033[0m"
        return 1
    fi
    if [ -s "$1" ]
    then
        echo -e "\033[32mPass\033[0m"
        return 0
    else
        echo -e "\033[31mFAILED. File is empty: $1\033[0m"
        return 1
    fi
}

assertEquals(){
    # USAGE
    # assertEquals $foo $bar
    exp=$1
    obs=$2
    if [ "$exp" == "$obs" ]
    then
        echo -e "\033[32mPass\033[0m"
        return 0
    else
        echo -e "\033[31mFAILED: expected '$exp' got '$obs'\033[0m"
    fi
    return 1
}

assertNotEquals(){
    # USAGE
    # assertNotEquals $foo $bar
    exp=$1
    obs=$2
    if [ "$exp" == "$obs" ]
    then
        echo -e "\033[31mFAILED: '$exp' is equal to '$obs'\033[0m"
        return 1
    else
        echo -e "\033[32mPass\033[0m"
        return 0
    fi
}

pprint(){
    # Pretty print a header message
    # USAGE
    # pprint Hello world
    echo -e "\033[1m" $@ "\033[0m"
}
