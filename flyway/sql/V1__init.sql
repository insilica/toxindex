-- V1__init.sql

CREATE TABLE user_groups (
    group_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    name VARCHAR(255) NOT NULL UNIQUE,
    description TEXT,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE users (
    user_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    email VARCHAR(255) NOT NULL UNIQUE,
    hashpw VARCHAR(255) NOT NULL,
    token VARCHAR(255),
    stripe_customer_id VARCHAR(255),
    email_verified BOOLEAN DEFAULT FALSE,
    group_id UUID REFERENCES user_groups(group_id) ON DELETE SET NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE user_links (
    link_id SERIAL PRIMARY KEY,
    user_id UUID REFERENCES users(user_id),
    link_token VARCHAR(255) NOT NULL,
    expiration TIMESTAMPTZ NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE workflows (
    workflow_id SERIAL PRIMARY KEY,
    title VARCHAR(255) NOT NULL UNIQUE,
    description TEXT,
    user_id UUID REFERENCES users(user_id) ON DELETE SET NULL,
    initial_prompt TEXT,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE workflow_group_access (
    access_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    workflow_id INTEGER REFERENCES workflows(workflow_id) ON DELETE CASCADE,
    group_id UUID REFERENCES user_groups(group_id) ON DELETE CASCADE,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(workflow_id, group_id)
);

CREATE TABLE environments (
    environment_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    title VARCHAR(255) NOT NULL,
    description TEXT,
    user_id UUID REFERENCES users(user_id) ON DELETE CASCADE,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE chat_sessions (
    session_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    environment_id UUID REFERENCES environments(environment_id) ON DELETE CASCADE,
    user_id UUID REFERENCES users(user_id) ON DELETE CASCADE,
    title TEXT,
    created_at TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE tasks (
    task_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    title VARCHAR(255) NOT NULL,
    description TEXT,
    user_id UUID REFERENCES users(user_id),
    workflow_id INTEGER REFERENCES workflows(workflow_id),
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    finished_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    modified_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    celery_task_id TEXT,
    environment_id UUID REFERENCES environments(environment_id),
    archived BOOLEAN DEFAULT FALSE,
    last_accessed TIMESTAMPTZ,
    session_id UUID REFERENCES chat_sessions(session_id) ON DELETE CASCADE,
    status VARCHAR(32) DEFAULT 'processing'
);

CREATE TABLE files (
    file_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    task_id UUID REFERENCES tasks(task_id) ON DELETE CASCADE,
    user_id UUID REFERENCES users(user_id),
    filename TEXT NOT NULL,
    filepath TEXT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    environment_id UUID REFERENCES environments(environment_id)
);

CREATE TABLE messages (
    message_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    task_id UUID REFERENCES tasks(task_id) ON DELETE CASCADE,
    user_id UUID REFERENCES users(user_id),
    role VARCHAR(50) NOT NULL,
    content TEXT NOT NULL,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    session_id UUID REFERENCES chat_sessions(session_id) ON DELETE CASCADE
);

-- Insert default user groups
INSERT INTO user_groups (group_id, name, description) VALUES
    ('00000000-0000-0000-0000-000000000001', 'admin', 'Administrators with full access'),
    ('00000000-0000-0000-0000-000000000002', 'researcher', 'Researchers with access to analysis workflows'),
    ('00000000-0000-0000-0000-000000000003', 'basic', 'Basic users with limited access')
ON CONFLICT (name) DO NOTHING;

-- Insert assistant user for system messages
INSERT INTO users (user_id, email, hashpw, token, stripe_customer_id)
VALUES (
  '00000000-0000-0000-0000-000000000000',
  'assistant@test.com',
  '', '', ''
)
ON CONFLICT (user_id) DO NOTHING;