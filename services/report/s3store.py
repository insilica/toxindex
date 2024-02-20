import boto3
import os
import dotenv
import logging
from botocore.exceptions import NoCredentialsError, ClientError

class S3Store:

    def __init__(self):
        dotenv.load_dotenv()
        self.access_key = os.environ.get('AWS_ACCESS_KEY')
        self.secret_key = os.environ.get('AWS_SECRET_ACCESS_KEY')
        self.bucket_name = os.environ.get('AWS_S3_BUCKET_NAME')
        self.s3_client = boto3.client('s3', aws_access_key_id=self.access_key, aws_secret_access_key=self.secret_key)

    def url_to_object(self, url):
        """Extract the object name from a URL

        :param url: URL to an S3 object
        :return: Object name
        """
        return url.split('/')[-1]
    
    def object_url(self, object_name):
        return f"https://{self.bucket_name}.s3.amazonaws.com/{object_name}"
    
    def upload_file(self, file_name, object_name=None):
        """Upload a file to an S3 bucket

        :param file_name: File to upload
        :param object_name: S3 object name. If not specified, file_name is used
        :return: True if file was uploaded, else False
        """
        if object_name is None:
            object_name = file_name

        try:
            self.s3_client.upload_file(file_name, self.bucket_name, object_name)
            logging.info(f"Uploaded {file_name} to {object_name}")
            return self.object_url(object_name)
        except ClientError as e:
            logging.error(e)
            return False
        except NoCredentialsError:
            logging.error("Credentials not available")
            return False

    def delete_object(self, object_name):
        """Delete an object from S3

        :param object_name: S3 object name
        :return: True if the referenced object was deleted, else False
        """
        try:
            self.s3_client.delete_object(Bucket=self.bucket_name, Key=object_name)
            logging.info(f"Deleted {object_name}")
            return True
        except ClientError as e:
            logging.error(e)
            return False
        except NoCredentialsError:
            logging.error("Credentials not available")
            return False

    def generate_presigned_url(self, object_name, expire=3600):
        """Generate a presigned URL to share an S3 object

        :param object_name: string
        :param expiration: Time in seconds for the presigned URL to remain valid
        :return: Presigned URL as string. If error, returns None.
        """
        try:
            p = {'Bucket': self.bucket_name, 'Key': object_name}
            response = self.s3_client.generate_presigned_url('get_object', Params=p, ExpiresIn=expire)
        except ClientError as e:
            logging.error(e)
            return None

        return response
    
    def download_file(self, object_name, file_name=None):
        """Download an object from S3

        :param object_name: S3 object name
        :param file_name: File name to save the object to. If not specified, object_name is used
        :return: True if file was downloaded, else False
        """
        if file_name is None:
            file_name = object_name

        try:
            self.s3_client.download_file(self.bucket_name, object_name, file_name)
            logging.info(f"Downloaded {object_name} to {file_name}")
            return file_name
        except ClientError as e:
            logging.error(e)
            return False
        except NoCredentialsError:
            logging.error("Credentials not available")
            return False
