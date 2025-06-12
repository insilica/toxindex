-- Add session_id to tasks table for chat session support
ALTER TABLE tasks ADD COLUMN session_id UUID REFERENCES chat_sessions(session_id);
-- Existing tasks will have session_id as NULL 