
# gives others access:
# make directories executable
find . -type d -exec chmod 755 {} \;

# make files readable
find . -type f -exec chmod 644 {} \;



# no "other" access:
 find . -type d -exec chmod 750 {} \; ; find . -type f -exec chmod 640 {} \;
