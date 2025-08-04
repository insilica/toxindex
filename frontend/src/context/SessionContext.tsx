import React, { createContext, useContext, useEffect, useRef, useState } from 'react';
import type { ReactNode } from 'react';

interface SessionSettings {
  timeoutMinutes: number;
  warningMinutes: number;
  refreshIntervalMinutes: number;
}

interface SessionStats {
  sessionStartTime: number;
  lastActivityTime: number;
  lastRefreshTime: number;
  refreshCount: number;
  totalSessionDuration: number;
  timeSinceLastActivity: number;
  timeSinceLastRefresh: number;
}

interface SessionState {
  // Session warning and timeout
  showWarning: boolean;
  timeRemaining: number;
  timeUntilTimeout: number;
  
  // Session settings
  settings: SessionSettings;
  settingsLoaded: boolean;
  
  // Session statistics
  stats: SessionStats;
  
  // Actions
  updateActivity: () => void;
  refreshSession: () => Promise<void>;
  extendSession: () => Promise<void>;
  logoutNow: () => void;
}

interface SessionProviderProps {
  children: ReactNode;
  onSessionExpired?: () => void;
}

const SessionContext = createContext<SessionState | undefined>(undefined);

export const SessionProvider: React.FC<SessionProviderProps> = ({ 
  children, 
  onSessionExpired 
}) => {
  const refreshIntervalRef = useRef<number | null>(null);
  const lastActivityRef = useRef<number>(Date.now());
  const lastRefreshRef = useRef<number>(Date.now());
  const sessionStartRef = useRef<number>(Date.now());
  const isRefreshingRef = useRef<boolean>(false);
  const onSessionExpiredRef = useRef(onSessionExpired);
  
  // Update the ref when onSessionExpired changes
  useEffect(() => {
    onSessionExpiredRef.current = onSessionExpired;
  }, [onSessionExpired]);
  
  const [showWarning, setShowWarning] = useState(false);
  const [timeRemaining, setTimeRemaining] = useState(0);
  const [refreshCount, setRefreshCount] = useState(0);
  
  const [settings, setSettings] = useState<SessionSettings>({
    timeoutMinutes: 60,
    warningMinutes: 5,
    refreshIntervalMinutes: 30
  });
  const [settingsLoaded, setSettingsLoaded] = useState(false);
  
  const [stats, setStats] = useState<SessionStats>({
    sessionStartTime: sessionStartRef.current,
    lastActivityTime: lastActivityRef.current,
    lastRefreshTime: lastRefreshRef.current,
    refreshCount: 0,
    totalSessionDuration: 0,
    timeSinceLastActivity: 0,
    timeSinceLastRefresh: 0
  });

  // Load settings from backend
  const loadSettings = async () => {
    try {
      console.log('[SessionContext] Loading session settings...');
      const response = await fetch('/api/auth/session_settings', {
        credentials: 'include',
      });

      console.log('[SessionContext] Settings response status:', response.status);
      
      if (response.ok) {
        const data = await response.json();
        console.log('[SessionContext] Settings response data:', data);
        if (data.success && data.settings) {
          console.log('[SessionContext] Loaded settings:', data.settings);
          setSettings({
            timeoutMinutes: data.settings.session_timeout_minutes || 15,
            warningMinutes: data.settings.session_warning_minutes || 14,
            refreshIntervalMinutes: data.settings.session_refresh_interval_minutes || 30
          });
          setSettingsLoaded(true);
        } else {
          console.log('[SessionContext] Invalid settings data, using defaults');
          setSettings({
            timeoutMinutes: 15,
            warningMinutes: 14,
            refreshIntervalMinutes: 30
          });
          setSettingsLoaded(true);
        }
      } else {
        console.log('[SessionContext] Settings endpoint failed, using defaults');
        setSettings({
          timeoutMinutes: 15,
          warningMinutes: 14,
          refreshIntervalMinutes: 30
        });
        setSettingsLoaded(true);
      }
    } catch (error) {
      console.error('[SessionContext] Error loading settings:', error);
      setSettings({
        timeoutMinutes: 15,
        warningMinutes: 14,
        refreshIntervalMinutes: 30
      });
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
    
    // Prevent concurrent refresh calls
    if (isRefreshingRef.current) {
      console.log('[SessionContext] Refresh already in progress, skipping');
      return;
    }
    
    if (timeSinceLastRefresh < minRefreshInterval) {
      console.log(`[SessionContext] Skipping refresh - too soon (${Math.round(timeSinceLastRefresh / 1000)}s since last refresh)`);
      return;
    }
    
    isRefreshingRef.current = true;
    
    try {
      console.log('[SessionContext] Sending refresh request...');
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
          console.log('[SessionContext] Session expired, logging out');
          onSessionExpiredRef.current?.();
        } else {
          console.log(`[SessionContext] Refresh failed with status: ${response.status}`);
        }
      } else {
        console.log('[SessionContext] Session refreshed successfully');
        lastRefreshRef.current = now;
        setRefreshCount(prev => prev + 1);
        setShowWarning(false);
      }
    } catch (error) {
      console.error('[SessionContext] Error refreshing session:', error);
    } finally {
      isRefreshingRef.current = false;
    }
  };

  // Extend session manually
  const extendSession = async () => {
    await refreshSession();
  };

  // Logout immediately
  const logoutNow = () => {
    onSessionExpiredRef.current?.();
  };

  // Calculate time until timeout
  const getTimeUntilTimeout = (): number => {
    if (!settings.timeoutMinutes || !stats.timeSinceLastActivity) {
      return 0;
    }
    const timeoutMs = settings.timeoutMinutes * 60 * 1000;
    return Math.max(0, timeoutMs - stats.timeSinceLastActivity);
  };

  // Update stats every second
  useEffect(() => {
    const interval = setInterval(() => {
      const now = Date.now();
      setStats(prev => ({
        ...prev,
        sessionStartTime: sessionStartRef.current,
        lastActivityTime: lastActivityRef.current,
        lastRefreshTime: lastRefreshRef.current,
        refreshCount,
        totalSessionDuration: now - sessionStartRef.current,
        timeSinceLastActivity: now - lastActivityRef.current,
        timeSinceLastRefresh: now - lastRefreshRef.current
      }));
    }, 1000);

    return () => clearInterval(interval);
  }, [refreshCount]);

  useEffect(() => {
    // Load settings on mount
    loadSettings();
  }, []);

  useEffect(() => {
    if (!settingsLoaded) return;

    // Ensure we have valid settings
    if (!settings.refreshIntervalMinutes || settings.refreshIntervalMinutes <= 0) {
      console.log('[SessionContext] Invalid refresh interval, using default of 30 minutes');
      setSettings(prev => ({ ...prev, refreshIntervalMinutes: 30 }));
      return;
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

      console.log(`[SessionContext] Checking session - Time since last activity: ${Math.round(timeSinceLastActivity / 1000)}s, Refresh interval: ${settings.refreshIntervalMinutes}min`);

      // If user has been inactive for too long, don't refresh
      if (timeSinceLastActivity > timeoutMs) {
        console.log('[SessionContext] User inactive for too long, not refreshing session');
        onSessionExpiredRef.current?.();
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
        console.log('[SessionContext] Refreshing session...');
        refreshSession();
      } else {
        console.log(`[SessionContext] Skipping refresh - ${Math.round((refreshIntervalMs - timeSinceLastRefresh) / 1000)}s until next refresh`);
      }
    }, 30000); // Check every 30 seconds

    refreshIntervalRef.current = refreshInterval;

    // Don't do initial refresh to prevent spam

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
  }, [settingsLoaded]);

  const contextValue: SessionState = {
    showWarning,
    timeRemaining,
    timeUntilTimeout: getTimeUntilTimeout(),
    settings,
    settingsLoaded,
    stats,
    updateActivity,
    refreshSession,
    extendSession,
    logoutNow,
  };

  return (
    <SessionContext.Provider value={contextValue}>
      {children}
    </SessionContext.Provider>
  );
};

export const useSession = (): SessionState => {
  const context = useContext(SessionContext);
  if (context === undefined) {
    throw new Error('useSession must be used within a SessionProvider');
  }
  return context;
};