import { useEffect, useRef, useState } from 'react';

interface UseSessionTimeoutOptions {
  onSessionExpired?: () => void;
}

export const useSessionTimeout = ({
  onSessionExpired
}: UseSessionTimeoutOptions = {}) => {
  const refreshIntervalRef = useRef<number | null>(null);
  const lastActivityRef = useRef<number>(Date.now());
  const lastRefreshRef = useRef<number>(0);
  const [showWarning, setShowWarning] = useState(false);
  const [timeRemaining, setTimeRemaining] = useState(0);
  const [settings, setSettings] = useState({
    timeoutMinutes: 60,
    warningMinutes: 5,
    refreshIntervalMinutes: 30
  });
  const [settingsLoaded, setSettingsLoaded] = useState(false);

  // Load settings from backend
  const loadSettings = async () => {
    try {
      const response = await fetch('/api/admin/settings/session', {
        credentials: 'include',
      });

      if (response.ok) {
        const data = await response.json();
        if (data.success) {
          console.log('[Session] Loaded settings:', data.settings);
          setSettings(data.settings);
          setSettingsLoaded(true);
        } else {
          console.log('[Session] Failed to load settings, using defaults');
          setSettingsLoaded(true);
        }
      } else {
        console.log('[Session] Settings endpoint failed, using defaults');
        setSettingsLoaded(true);
      }
    } catch (error) {
      console.error('[Session] Error loading settings:', error);
      // Use defaults if settings can't be loaded
      setSettingsLoaded(true);
    }
  };

  // Update last activity on user interaction
  const updateActivity = () => {
    lastActivityRef.current = Date.now();
    setShowWarning(false); // Hide warning when user is active
  };

  // Refresh session
  const refreshSession = async () => {
    const now = Date.now();
    const timeSinceLastRefresh = now - lastRefreshRef.current;
    const minRefreshInterval = 5000; // Minimum 5 seconds between refreshes
    
    if (timeSinceLastRefresh < minRefreshInterval) {
      console.log(`[Session] Skipping refresh - too soon (${Math.round(timeSinceLastRefresh / 1000)}s since last refresh)`);
      return;
    }
    
    try {
      console.log('[Session] Sending refresh request...');
      const response = await fetch('/api/auth/refresh_session', {
        method: 'POST',
        credentials: 'include',
        headers: {
          'Content-Type': 'application/json',
        },
      });

      if (!response.ok) {
        if (response.status === 401) {
          // Session expired
          console.log('[Session] Session expired, logging out');
          onSessionExpired?.();
        } else {
          console.log(`[Session] Refresh failed with status: ${response.status}`);
        }
      } else {
        console.log('[Session] Session refreshed successfully');
        lastRefreshRef.current = now;
        setShowWarning(false);
      }
    } catch (error) {
      console.error('[Session] Error refreshing session:', error);
    }
  };

  // Extend session manually
  const extendSession = async () => {
    await refreshSession();
  };

  // Logout immediately
  const logoutNow = () => {
    onSessionExpired?.();
  };

  useEffect(() => {
    // Load settings on mount
    loadSettings();
  }, []);

  useEffect(() => {
    if (!settingsLoaded) return;

    // Ensure we have valid settings
    if (!settings.refreshIntervalMinutes || settings.refreshIntervalMinutes <= 0) {
      console.log('[Session] Invalid refresh interval, using default of 30 minutes');
      setSettings(prev => ({ ...prev, refreshIntervalMinutes: 30 }));
      return; // This will trigger the effect again with valid settings
    }

    // Clean up any existing interval first
    if (refreshIntervalRef.current) {
      clearInterval(refreshIntervalRef.current);
      refreshIntervalRef.current = null;
    }

    // Set up activity listeners
    const events = ['mousedown', 'mousemove', 'keypress', 'scroll', 'touchstart', 'click'];
    
    const handleActivity = () => {
      updateActivity();
    };

    events.forEach(event => {
      document.addEventListener(event, handleActivity, true);
    });

    // Set up periodic session refresh and warning checks
    const refreshInterval = setInterval(() => {
      const now = Date.now();
      const timeSinceLastActivity = now - lastActivityRef.current;
      const timeoutMs = settings.timeoutMinutes * 60 * 1000;
      const warningMs = settings.warningMinutes * 60 * 1000;
      const refreshIntervalMs = settings.refreshIntervalMinutes * 60 * 1000;

      console.log(`[Session] Checking session - Time since last activity: ${Math.round(timeSinceLastActivity / 1000)}s, Refresh interval: ${settings.refreshIntervalMinutes}min`);

      // If user has been inactive for too long, don't refresh
      if (timeSinceLastActivity > timeoutMs) {
        console.log('[Session] User inactive for too long, not refreshing session');
        onSessionExpired?.();
        return;
      }

      // Check if we should show warning
      const timeUntilExpiry = timeoutMs - timeSinceLastActivity;
      if (timeUntilExpiry <= warningMs && timeUntilExpiry > 0) {
        setShowWarning(true);
        setTimeRemaining(Math.ceil(timeUntilExpiry / 1000));
      } else {
        setShowWarning(false);
      }

      // Only refresh if enough time has passed since last refresh
      const timeSinceLastRefresh = now - lastRefreshRef.current;
      if (timeSinceLastRefresh >= refreshIntervalMs) {
        console.log('[Session] Refreshing session...');
        refreshSession();
      } else {
        console.log(`[Session] Skipping refresh - ${Math.round((refreshIntervalMs - timeSinceLastRefresh) / 1000)}s until next refresh`);
      }
    }, Math.min(60000, settings.refreshIntervalMinutes * 60 * 1000)); // Check every minute or refresh interval, whichever is smaller

    refreshIntervalRef.current = refreshInterval;

    // Initial session refresh
    refreshSession();

    return () => {
      // Clean up event listeners
      events.forEach(event => {
        document.removeEventListener(event, handleActivity, true);
      });

      // Clean up interval
      if (refreshIntervalRef.current) {
        clearInterval(refreshIntervalRef.current);
        refreshIntervalRef.current = null;
      }
    };
  }, [settingsLoaded, onSessionExpired]); // Removed 'settings' from dependency array

  return {
    updateActivity,
    refreshSession,
    extendSession,
    logoutNow,
    showWarning,
    timeRemaining,
  };
}; 