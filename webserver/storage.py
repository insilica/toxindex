import os
from typing import Optional
import boto3


class S3FileStorage:
    """Simple wrapper around boto3 S3 client for uploading and downloading files."""

    def __init__(self, bucket_name: Optional[str] = None, region_name: Optional[str] = None):
        self.bucket_name = bucket_name or os.getenv("S3_BUCKET")
        self.region_name = region_name or os.getenv("AWS_REGION")
        if not self.bucket_name:
            raise ValueError("S3 bucket name must be provided via argument or S3_BUCKET env var")
        self.s3 = boto3.client("s3", region_name=self.region_name)

    def upload_file(self, file_path: str, key: Optional[str] = None, content_type: Optional[str] = None) -> str:
        """Upload a local file to S3 and return the object key."""
        if key is None:
            key = os.path.basename(file_path)
        extra_args = {"ContentType": content_type} if content_type else None
        if extra_args:
            self.s3.upload_file(file_path, self.bucket_name, key, ExtraArgs=extra_args)
        else:
            self.s3.upload_file(file_path, self.bucket_name, key)
        return key

    def generate_download_url(self, key: str, expires_in: int = 3600) -> str:
        """Generate a presigned URL to download a file."""
        return self.s3.generate_presigned_url(
            ClientMethod="get_object",
            Params={"Bucket": self.bucket_name, "Key": key},
            ExpiresIn=expires_in,
        )
