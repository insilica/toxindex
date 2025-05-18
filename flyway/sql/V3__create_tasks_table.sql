CREATE TABLE tasks (
    task_id SERIAL PRIMARY KEY,
    title VARCHAR(255) NOT NULL,
    description TEXT,
    user_id INTEGER REFERENCES users(user_id),
    workflow_id INTEGER REFERENCES workflows(workflow_id),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);