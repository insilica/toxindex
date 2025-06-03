ALTER TABLE messages DROP CONSTRAINT messages_task_id_fkey;
ALTER TABLE messages
ADD CONSTRAINT messages_task_id_fkey
FOREIGN KEY (task_id) REFERENCES tasks(task_id) ON DELETE CASCADE; 