-- UUID migration for all IDs and foreign keys (was V14, now V13)

-- 0. Ensure chat_sessions table and session_id column exist
CREATE TABLE IF NOT EXISTS chat_sessions (
    session_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    environment_id INTEGER REFERENCES environments(environment_id) ON DELETE CASCADE,
    user_id UUID REFERENCES users(user_id) ON DELETE CASCADE,
    title TEXT,
    created_at TIMESTAMP DEFAULT now()
);

DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_name='messages' AND column_name='session_id'
    ) THEN
        ALTER TABLE messages ADD COLUMN session_id UUID REFERENCES chat_sessions(session_id);
    END IF;
END$$;

-- 1. Add new UUID columns
ALTER TABLE environments ADD COLUMN environment_id_new UUID DEFAULT gen_random_uuid();
ALTER TABLE tasks ADD COLUMN task_id_new UUID DEFAULT gen_random_uuid();
ALTER TABLE tasks ADD COLUMN environment_id_new UUID;
ALTER TABLE files ADD COLUMN file_id_new UUID DEFAULT gen_random_uuid();
ALTER TABLE files ADD COLUMN environment_id_new UUID;
ALTER TABLE files ADD COLUMN task_id_new UUID;
ALTER TABLE chat_sessions ADD COLUMN session_id_new UUID DEFAULT gen_random_uuid();
ALTER TABLE chat_sessions ADD COLUMN environment_id_new UUID;
ALTER TABLE messages ADD COLUMN message_id_new UUID DEFAULT gen_random_uuid();
ALTER TABLE messages ADD COLUMN task_id_new UUID;
ALTER TABLE messages ADD COLUMN session_id_new UUID;

-- 2. Populate UUIDs for all existing rows
UPDATE environments SET environment_id_new = gen_random_uuid();
UPDATE tasks SET task_id_new = gen_random_uuid();
UPDATE files SET file_id_new = gen_random_uuid();
UPDATE chat_sessions SET session_id_new = gen_random_uuid();
UPDATE messages SET message_id_new = gen_random_uuid();

-- 3. Update foreign key references
-- tasks.environment_id_new
UPDATE tasks SET environment_id_new = e.environment_id_new FROM environments e WHERE tasks.environment_id = e.environment_id;
-- files.environment_id_new, files.task_id_new
UPDATE files SET environment_id_new = e.environment_id_new FROM environments e WHERE files.environment_id = e.environment_id;
UPDATE files SET task_id_new = t.task_id_new FROM tasks t WHERE files.task_id = t.task_id;
-- chat_sessions.environment_id_new
UPDATE chat_sessions SET environment_id_new = e.environment_id_new FROM environments e WHERE chat_sessions.environment_id = e.environment_id;
-- messages.task_id_new, messages.session_id_new
UPDATE messages SET task_id_new = t.task_id_new FROM tasks t WHERE messages.task_id = t.task_id;
UPDATE messages SET session_id_new = s.session_id_new FROM chat_sessions s WHERE messages.session_id = s.session_id;

-- 4. Drop old constraints (if any)
ALTER TABLE tasks DROP CONSTRAINT IF EXISTS tasks_environment_id_fkey;
ALTER TABLE files DROP CONSTRAINT IF EXISTS files_environment_id_fkey;
ALTER TABLE files DROP CONSTRAINT IF EXISTS files_task_id_fkey;
ALTER TABLE chat_sessions DROP CONSTRAINT IF EXISTS chat_sessions_environment_id_fkey;
ALTER TABLE messages DROP CONSTRAINT IF EXISTS messages_task_id_fkey;
ALTER TABLE messages DROP CONSTRAINT IF EXISTS messages_session_id_fkey;

-- 5. Drop old integer columns and rename new columns
ALTER TABLE environments DROP COLUMN environment_id;
ALTER TABLE environments RENAME COLUMN environment_id_new TO environment_id;

-- Drop the primary key constraint by name if it exists
DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM information_schema.table_constraints
        WHERE table_name = 'tasks' AND constraint_type = 'PRIMARY KEY'
    ) THEN
        ALTER TABLE tasks DROP CONSTRAINT tasks_pkey;
    END IF;
END$$;

ALTER TABLE tasks DROP COLUMN environment_id;
ALTER TABLE tasks DROP COLUMN task_id;
ALTER TABLE tasks RENAME COLUMN environment_id_new TO environment_id;
ALTER TABLE tasks RENAME COLUMN task_id_new TO task_id;
ALTER TABLE tasks ADD PRIMARY KEY (task_id);
ALTER TABLE files DROP COLUMN environment_id;
ALTER TABLE files DROP COLUMN task_id;
ALTER TABLE files DROP COLUMN file_id;
ALTER TABLE files RENAME COLUMN environment_id_new TO environment_id;
ALTER TABLE files RENAME COLUMN task_id_new TO task_id;
ALTER TABLE files RENAME COLUMN file_id_new TO file_id;
ALTER TABLE chat_sessions DROP CONSTRAINT IF EXISTS chat_sessions_pkey;
ALTER TABLE chat_sessions DROP COLUMN environment_id;
ALTER TABLE chat_sessions DROP COLUMN session_id;
ALTER TABLE chat_sessions RENAME COLUMN environment_id_new TO environment_id;
ALTER TABLE chat_sessions RENAME COLUMN session_id_new TO session_id;
ALTER TABLE chat_sessions ADD PRIMARY KEY (session_id);
ALTER TABLE messages DROP COLUMN task_id;
ALTER TABLE messages DROP COLUMN session_id;
ALTER TABLE messages DROP COLUMN message_id;
ALTER TABLE messages RENAME COLUMN task_id_new TO task_id;
ALTER TABLE messages RENAME COLUMN session_id_new TO session_id;
ALTER TABLE messages RENAME COLUMN message_id_new TO message_id;

-- 6. Add new constraints
ALTER TABLE environments ADD PRIMARY KEY (environment_id);
ALTER TABLE tasks ADD PRIMARY KEY (task_id);
ALTER TABLE chat_sessions ADD PRIMARY KEY (session_id);
ALTER TABLE tasks ADD CONSTRAINT tasks_environment_id_fkey FOREIGN KEY (environment_id) REFERENCES environments(environment_id) ON DELETE CASCADE;
ALTER TABLE files ADD CONSTRAINT files_environment_id_fkey FOREIGN KEY (environment_id) REFERENCES environments(environment_id) ON DELETE CASCADE;
ALTER TABLE files ADD CONSTRAINT files_task_id_fkey FOREIGN KEY (task_id) REFERENCES tasks(task_id) ON DELETE CASCADE;
ALTER TABLE chat_sessions ADD CONSTRAINT chat_sessions_environment_id_fkey FOREIGN KEY (environment_id) REFERENCES environments(environment_id) ON DELETE CASCADE;
ALTER TABLE messages ADD CONSTRAINT messages_task_id_fkey FOREIGN KEY (task_id) REFERENCES tasks(task_id) ON DELETE CASCADE;
ALTER TABLE messages ADD CONSTRAINT messages_session_id_fkey FOREIGN KEY (session_id) REFERENCES chat_sessions(session_id) ON DELETE CASCADE;

-- Add session_id to tasks table for chat session support
ALTER TABLE tasks ADD COLUMN session_id UUID REFERENCES chat_sessions(session_id); 