ssh jon@ki-rh42.slac.stanford.edu -T -c arcfour -o Compression=no -x "tar cf - /u1/jon" | tar xf - -C .
