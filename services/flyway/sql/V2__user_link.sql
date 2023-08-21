CREATE TABLE user_link (
 user_link_id SERIAL PRIMARY KEY,
 link_token TEXT NOT NULL UNIQUE,
 expiration TIMESTAMP NOT NULL,
 user_id INTEGER,
 FOREIGN KEY (user_id) REFERENCES users(user_id) ON DELETE CASCADE
);
