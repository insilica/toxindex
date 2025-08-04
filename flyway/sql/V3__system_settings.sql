-- V3__system_settings.sql

CREATE TABLE system_settings (
    setting_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    setting_key VARCHAR(255) NOT NULL UNIQUE,
    setting_value TEXT NOT NULL,
    description TEXT,
    created_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP
);

-- Insert default system settings
INSERT INTO system_settings (setting_key, setting_value, description) VALUES
    ('session_timeout_minutes', '60', 'Session timeout in minutes for user authentication'),
    ('session_warning_minutes', '5', 'Warning time in minutes before session expires'),
    ('session_refresh_interval_minutes', '30', 'Interval in minutes for automatic session refresh')
ON CONFLICT (setting_key) DO NOTHING;

-- Create index for faster lookups
CREATE INDEX idx_system_settings_key ON system_settings(setting_key); 