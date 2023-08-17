CREATE TABLE users (
 user_id   INTEGER PRIMARY KEY,
 email TEXT NOT NULL UNIQUE,
 token TEXT NOT NULL UNIQUE,
 stripe_customer_id TEXT NOT NULL UNIQUE,
 hashpw TEXT NOT NULL,
 email_verified BOOL DEFAULT FALSE CHECK (email_verified IN (FALSE, TRUE)),
 registration_date TIMESTAMP DEFAULT current_timestamp
);
