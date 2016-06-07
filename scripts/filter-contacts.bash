#!/bin/bash

if [ ! $2 ] ; then
    echo "Usage: ./script <chain groups to ignore> <contacts file>"
    echo "       example: ./remove-chain-contacts.bash \"ABC,DE\" test01.contacts"
    echo "       The above example removes contacts between chains A, B and C, as well as betwen D and E."
    exit
fi

str=$1
ifile=$2

# CentOS awk/gawk does NOT come with chr() for some reason.
awk -v igChainStr=$str \
'function chToChar(c)
{
    return sprintf("%c", c + 64)
}

BEGIN {
    #Interpret chain string
    split(igChainStr, ig_txtarr, ",")
    #for (i in ig_txtarr) {print ig_txtarr[i]}
}
{
    segA=chToChar(int(substr($0,1,2)))
    segB=chToChar(int(substr($0,3,2)))
    #print "= = Contact between chain " segA " and " segB
    for (group in ig_txtarr) {
        #print group, segA, match(ig_txtarr[group], segA), segB, match(ig_txtarr[group], segB)
        if ( match(ig_txtarr[group], segA) > 0 && match(ig_txtarr[group], segB) > 0 ) {
            #print "= = Ignoring this pair as it is found in symmetry group " group
            next
        }
    }
    print $0
}' $ifile
