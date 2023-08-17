CREATE TABLE user_link(
 ulink Integer PRIMARY KEY,
 uid   INTEGER NOT NULL,
 link_token TEXT NOT NULL UNIQUE,
 expiration DATE NOT NULL
)