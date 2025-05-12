CREATE TABLE messages (
    message_id SERIAL PRIMARY KEY,
    run_id INTEGER REFERENCES runs(run_id),
    user_id INTEGER REFERENCES users(user_id),
    role VARCHAR(50) NOT NULL,
    content TEXT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);