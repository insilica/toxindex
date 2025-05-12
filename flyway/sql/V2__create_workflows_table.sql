CREATE TABLE workflows (
    workflow_id SERIAL PRIMARY KEY,
    title VARCHAR(255) NOT NULL UNIQUE,
    description TEXT,
    user_id INTEGER REFERENCES users(user_id) NULL,

    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);