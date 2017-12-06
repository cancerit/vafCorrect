echo "about to run $( which perl ) which is $( readlink -f "$( which perl )")"
echo "about to run $( which prove ) which is $( readlink -f "$( which prove )")"
prove -v -I lib 
