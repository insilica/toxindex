CREATE TABLE files (
    file_id SERIAL PRIMARY KEY,
    task_id INTEGER REFERENCES tasks(task_id),
    user_id UUID REFERENCES users(user_id),
    filename TEXT NOT NULL,
    filepath TEXT NOT NULL,        -- path within the bucket or filesystem
    s3_url TEXT NOT NULL,          -- full public or signed S3 URL
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
