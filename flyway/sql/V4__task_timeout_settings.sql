-- V4__task_timeout_settings.sql

-- Add task timeout settings to system_settings
INSERT INTO system_settings (setting_key, setting_value, description) VALUES
    ('task_timeout_minutes', '{"toxindex_rap": 10, "toxindex_vanilla": 10, "toxindex_json": 10, "raptool": 10, "pathway_analysis": 10, "default": 10}', 'Task timeout settings in minutes for different workflows'),
    ('task_timeout_toxindex_rap', '10', 'Timeout in minutes for ToxIndex RAP workflow'),
    ('task_timeout_toxindex_vanilla', '10', 'Timeout in minutes for ToxIndex Vanilla workflow'),
    ('task_timeout_toxindex_json', '10', 'Timeout in minutes for ToxIndex JSON workflow'),
    ('task_timeout_raptool', '10', 'Timeout in minutes for RAPtool workflow'),
    ('task_timeout_pathway_analysis', '10', 'Timeout in minutes for Pathway Analysis workflow'),
    ('task_timeout_default', '10', 'Default timeout in minutes for workflows')
ON CONFLICT (setting_key) DO NOTHING;
