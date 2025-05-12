scp services/prod.env toxindex:/mnt/ebs/toxindex/services
scp services/report/data/cvae.sqlite toxindex:/mnt/ebs/toxindex/services/report/data/cvae.sqlite

# zip the entire directory
zip -r deploy.zip .

# scp the zip file to the server
scp deploy.zip toxindex:/home/ubuntu

# ssh into the server
ssh toxindex

# stop all running containers
docker stop $(docker ps -a -q)

# mv and unzip deploy.zip
mv ~/deploy.zip /mnt/ebs/toxindex
cd /mnt/ebs


# unzip the file to /mnt/toxindex
unzip deploy.zip 