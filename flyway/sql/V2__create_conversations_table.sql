CREATE TABLE conversations (
    conversation_id SERIAL PRIMARY KEY,
    title VARCHAR(255) NOT NULL,
    user_id INTEGER REFERENCES users(user_id),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);