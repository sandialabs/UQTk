#!/bin/bash
read -p "enter directory: " directory

echo "directory = $directory"
#file that all changed names get saved to
FILECHANGED="files_changed.new"
FILE_NOT_CHANGED="files_not_changed.new"

#for all the python files and shell scripts
for i in $directory/*.py $directory/*.sh; do
    #get header file name and check that it is there
    LICENCEFILE="new_header_py"
    [ ! -f "$LICENCEFILE" ] && echo "$LICENCEFILE is missing. Abort." && exit 1

    #If there is no header
    if [ "$(head -c2 $i)" == "#=" ]
    then
      echo "$i being converted"
      NEWFILE="${i}.new"
      [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

      cat "$LICENCEFILE" >> "$NEWFILE"
      sed '1,26d' "$i" >> "$NEWFILE"
      echo "$i" >> $FILECHANGED

    #if needed copy the first line then do header
  elif [ "$(head -2 $i |tail -1 | head -c2)" == "#=" ]
    then
        echo "$i being converted with a 1 line header"
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

        head -1 $i >> "$NEWFILE"
        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,27d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

    elif [ "$(head -3 $i |tail -2 | head -c2)" == "#=" ]
    then
        echo "$i being converted with a 2 line header"
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

        head -2 $i >> "$NEWFILE"
        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,28d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

    else
      echo "Could not find start (#=) for file $i."
      echo "$i" >> $FILE_NOT_CHANGED

    fi

done

#for all the .i files

for i in $directory/*.i; do
    #get header file name and check that it is there
    LICENCEFILE="new_header_i"
    [ ! -f "$LICENCEFILE" ] && echo "$LICENCEFILE is missing. Abort." && exit 1

    #If there is no header
    if [ "$(head -c3 $i)" == "//=" ]
    then
      echo "$i being converted"
      NEWFILE="${i}.new"
      [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

      cat "$LICENCEFILE" >> "$NEWFILE"
      sed '1,26d' "$i" >> "$NEWFILE"
      echo "$i" >> $FILECHANGED

    #if needed copy the first line then do header
  elif [ "$(head -2 $i |tail -1 | head -c3)" == "//=" ]
    then
        echo "$i being converted with a 1 line header"
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

        head -1 $i >> "$NEWFILE"
        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,27d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

    elif [ "$(head -3 $i |tail -2 | head -c3)" == "//=" ]
    then
        echo "$i being converted with a 2 line header"
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

        head -2 $i >> "$NEWFILE"
        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,28d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

    else
      echo "Could not find start (//=) for file $i."
      echo "$i" >> $FILE_NOT_CHANGED

    fi

done


#now check all .cpp file and .h files
for i in $directory/*.cpp $directory/*.h; do
    #get header file name and check that it is there
    LICENCEFILE="new_header_cpp"
    [ ! -f "$LICENCEFILE" ] && echo "$LICENCEFILE is missing. Abort." && exit 1

    #If there is no header
    if [ "$(head -c4 $i)" == "/* =" ]
    then
      echo "$i being converted"
      NEWFILE="${i}.new"
      [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

      cat "$LICENCEFILE" >> "$NEWFILE"
      sed '1,26d' "$i" >> "$NEWFILE"
      echo "$i" >> $FILECHANGED

    #if needed copy the first line then do header
  elif [ "$(head -2 $i |tail -1 | head -c4)" == "/* ="  ]
    then
        echo "$i being converted with a 1 line header"
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

        head -1 $i >> "$NEWFILE"
        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,27d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

    elif [ "$(head -3 $i |tail -2 | head -c4)" == "/* ="  ]
    then
        echo "$i being converted with a 2 line header"
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue

        head -2 $i >> "$NEWFILE"
        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,28d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

    else
      echo "Could not find start (/* =) for file $i."
      echo "$i" >> $FILE_NOT_CHANGED
    fi

done
