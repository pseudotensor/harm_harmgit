sudo bash
name="career"
mkdir /home/svn/$name
svnadmin create /home/svn/$name
cd /home/svn/$name/
mkdir trunk tags branches
chown -R www-data:www-data .
chmod -R g+rws .
cd /data/jon/
svn import $name file:///home/svn/$name -m 'initial import'
cd /home/svn/$name/
chown -R www-data:www-data .
chmod -R g+rws .
exit

