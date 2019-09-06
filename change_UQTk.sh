#!/bin/bash

#enter the name of the directory, so we can search through all directories manually pretty easily
read -p "enter directory: " directory_list

echo "Changing from directory: $directory_list"
#file that all changed names get saved to
FILECHANGED="files_changed.new"
FILE_NOT_CHANGED="files_not_changed.new"

for directory in $directory_list; do
  #for all the python files and shell scripts
  for i in $directory/*.py $directory/*.sh; do
      #get header file name and check that it is there
      LICENCEFILE="new_header_py"
      [ ! -f "$LICENCEFILE" ] && echo "$LICENCEFILE is missing. Abort." && exit 1

      #If there is no header
      if [ "$(head -c2 $i)" == "#=" ] || [ "$(head -c3 $i)" == "# =" ]
      then
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
        echo "$i being converted"

        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,26d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

      #if needed copy the first line then do header
    elif [ "$(head -2 $i |tail -1 | head -c2)" == "#=" ] || [ "$(head -2 $i |tail -1 | head -c3)" == "# =" ]
      then
          NEWFILE="${i}.new"
          [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
          echo "$i being converted with a 1 line header"

          head -1 $i >> "$NEWFILE"
          cat "$LICENCEFILE" >> "$NEWFILE"
          sed '1,27d' "$i" >> "$NEWFILE"
          echo "$i" >> $FILECHANGED

      #if the header is 2 lines.
      elif [ "$(head -3 $i |tail -2 | head -c2)" == "#=" ]
      then
          NEWFILE="${i}.new"
          [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
          echo "$i being converted with a 2 line header"

          head -2 $i >> "$NEWFILE"
          cat "$LICENCEFILE" >> "$NEWFILE"
          sed '1,28d' "$i" >> "$NEWFILE"
          echo "$i" >> $FILECHANGED

      else
        echo "Could not find start (#=) or (# =)for file $i."
        if test -f "$i"; then
          echo "$i" >> $FILE_NOT_CHANGED
        fi

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
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
        echo "$i being converted"

        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,26d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

      #if needed copy the first line then do header
    elif [ "$(head -2 $i |tail -1 | head -c3)" == "//=" ]
      then
          NEWFILE="${i}.new"
          [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
          echo "$i being converted with a 1 line header"

          head -1 $i >> "$NEWFILE"
          cat "$LICENCEFILE" >> "$NEWFILE"
          sed '1,27d' "$i" >> "$NEWFILE"
          echo "$i" >> $FILECHANGED

      elif [ "$(head -3 $i |tail -2 | head -c3)" == "//=" ]
      then
          NEWFILE="${i}.new"
          [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
          echo "$i being converted with a 2 line header"

          head -2 $i >> "$NEWFILE"
          cat "$LICENCEFILE" >> "$NEWFILE"
          sed '1,28d' "$i" >> "$NEWFILE"
          echo "$i" >> $FILECHANGED

      else
        echo "Could not find start (//=) for file $i."
        if test -f "$i"; then
          echo "$i" >> $FILE_NOT_CHANGED
        fi

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
        NEWFILE="${i}.new"
        [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
        echo "$i being converted"

        cat "$LICENCEFILE" >> "$NEWFILE"
        sed '1,26d' "$i" >> "$NEWFILE"
        echo "$i" >> $FILECHANGED

      #if needed copy the first line then do header
    elif [ "$(head -2 $i |tail -1 | head -c4)" == "/* ="  ]
      then
          NEWFILE="${i}.new"
          [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
          echo "$i being converted with a 1 line header"

          head -1 $i >> "$NEWFILE"
          cat "$LICENCEFILE" >> "$NEWFILE"
          sed '1,27d' "$i" >> "$NEWFILE"
          echo "$i" >> $FILECHANGED

      elif [ "$(head -3 $i |tail -2 | head -c4)" == "/* ="  ]
      then
          NEWFILE="${i}.new"
          [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
          echo "$i being converted with a 2 line header"

          head -2 $i >> "$NEWFILE"
          cat "$LICENCEFILE" >> "$NEWFILE"
          sed '1,28d' "$i" >> "$NEWFILE"
          echo "$i" >> $FILECHANGED

      else
        echo "Could not find start (/* =) for file $i."
        if test -f "$i"; then
          echo "$i" >> $FILE_NOT_CHANGED
        fi
      fi
    done

    #now check all .f files (there is at least one in app/gkpSparse )
    for i in $directory/*.f; do
        #get header file name and check that it is there
        LICENCEFILE="new_header_f"
        [ ! -f "$LICENCEFILE" ] && echo "$LICENCEFILE is missing. Abort." && exit 1

        #If there is no header
        if [ "$(head -c5 $i)" == "c$$$=" ]
        then
          NEWFILE="${i}.new"
          [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
          echo "$i being converted"

          cat "$LICENCEFILE" >> "$NEWFILE"
          sed '1,26d' "$i" >> "$NEWFILE"
          echo "$i" >> $FILECHANGED

        #if needed copy the first line then do header
      elif [ "$(head -2 $i |tail -1 | head -c5)" == "c$$$="  ]
        then
            NEWFILE="${i}.new"
            [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
            echo "$i being converted with a 1 line header"

            head -1 $i >> "$NEWFILE"
            cat "$LICENCEFILE" >> "$NEWFILE"
            sed '1,27d' "$i" >> "$NEWFILE"
            echo "$i" >> $FILECHANGED

        elif [ "$(head -3 $i |tail -2 | head -c5)" == "c$$$="  ]
        then
            NEWFILE="${i}.new"
            [ -f "$NEWFILE" ] && echo "Sorry, $NEWFILE already exists" && continue
            echo "$i being converted with a 2 line header"

            head -2 $i >> "$NEWFILE"
            cat "$LICENCEFILE" >> "$NEWFILE"
            sed '1,28d' "$i" >> "$NEWFILE"
            echo "$i" >> $FILECHANGED

        else
          echo "Could not find start (c$$$=) for file $i."
          if test -f "$i"; then
            echo "$i" >> $FILE_NOT_CHANGED
          fi
        fi


      done
done
