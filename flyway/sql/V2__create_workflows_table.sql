CREATE TABLE workflows (
    workflow_id SERIAL PRIMARY KEY,
    title VARCHAR(255) NOT NULL UNIQUE,
    description TEXT,
    user_id UUID REFERENCES users(user_id) ON DELETE SET NULL,
    initial_prompt TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);