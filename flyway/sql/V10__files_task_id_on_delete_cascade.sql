-- Drop the existing foreign key constraint if it exists
ALTER TABLE files
DROP CONSTRAINT IF EXISTS files_task_id_fkey;

-- Add the foreign key constraint with ON DELETE CASCADE
ALTER TABLE files
ADD CONSTRAINT files_task_id_fkey
FOREIGN KEY (task_id) REFERENCES tasks(task_id) ON DELETE CASCADE; 