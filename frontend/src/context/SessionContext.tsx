import React, { createContext, useContext, useEffect, useRef, useState, useMemo, useCallback } from 'react';
import type { ReactNode } from 'react';

interface SessionSettings {
  timeoutMinutes: number;
  warningMinutes: number;
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

// Server-side session management component
const SessionManager: React.FC<{ children: ReactNode }> = ({ children }) => {
  const [showWarning, setShowWarning] = useState(false);
  const [timeRemaining, setTimeRemaining] = useState(0);
  const [isExpired, setIsExpired] = useState(false);
  const [sessionSettings, setSessionSettings] = useState({
    timeoutMinutes: 15,
    warningMinutes: 14
  });
  const statusCheckIntervalRef = useRef<number | null>(null);
  const activityUpdateIntervalRef = useRef<number | null>(null);

  const lastServerTimeRemainingRef = useRef<number>(0);
  const lastActivityUpdateRef = useRef<number>(0);
  const sessionSettingsRef = useRef(sessionSettings);

  // Update ref when settings change
  useEffect(() => {
    sessionSettingsRef.current = sessionSettings;
  }, [sessionSettings]);

  // Load session settings from server
  const loadSessionSettings = useCallback(async () => {
    try {
      const response = await fetch('/api/auth/session_settings', {
        credentials: 'include',
      });

      if (response.ok) {
        const data = await response.json();
        if (data.success && data.settings) {
          setSessionSettings({
            timeoutMinutes: data.settings.session_timeout_minutes || 15,
            warningMinutes: data.settings.session_warning_minutes || 14
          });
        }
      }
    } catch (error) {
      console.error('[SessionManager] Error loading session settings:', error);
    }
  }, []);

  // Check session status from server
  const checkSessionStatus = useCallback(async () => {
    try {
      const response = await fetch('/api/auth/session_status', {
        credentials: 'include',
      });

      if (!response.ok) {
        if (response.status === 401) {
          console.log('[SessionManager] Session expired on server');
          setIsExpired(true);
          setShowWarning(false);
          setTimeRemaining(0);
          return;
        }
        return;
      }

      const data = await response.json();

      if (data.success) {
        // Only update if values actually changed
        if (data.timeRemaining !== lastServerTimeRemainingRef.current) {
          lastServerTimeRemainingRef.current = data.timeRemaining;
          setTimeRemaining(data.timeRemaining);
        }
        
        if (data.showWarning !== showWarning) {
          setShowWarning(data.showWarning);
        }
        
        if (data.isExpired !== isExpired) {
          setIsExpired(data.isExpired);
        }
        
      }
    } catch (error) {
      console.error('[SessionManager] Error checking session status:', error);
    }
  }, [showWarning, isExpired]);

  // Update activity on server with dynamic interval based on session settings
  const updateActivity = useCallback(async () => {
    const now = Date.now();
    const timeSinceLastUpdate = now - lastActivityUpdateRef.current;
    
    // Calculate dynamic interval: session timeout minus warning time
    // This ensures we refresh before the warning appears
    const dynamicInterval = (sessionSettingsRef.current.timeoutMinutes - sessionSettingsRef.current.warningMinutes) * 60 * 1000;
    const minUpdateInterval = Math.max(dynamicInterval, 30000); // Minimum 30 seconds
    
    if (timeSinceLastUpdate < minUpdateInterval) {
      return; // Skip if too soon
    }
    
    try {
      lastActivityUpdateRef.current = now;
      console.log('[SessionManager] Updating activity on server...');
      
      await fetch('/api/auth/update_activity', {
        method: 'POST',
        credentials: 'include',
      });
      
      // Immediately check session status after updating activity
      // This ensures the countdown reflects the activity right away
      await checkSessionStatus();
    } catch (error) {
      console.error('[SessionManager] Error updating activity:', error);
    }
  }, [checkSessionStatus]);

  // Refresh session
  const refreshSession = useCallback(async () => {
    try {
      const response = await fetch('/api/auth/refresh_session', {
        method: 'POST',
        credentials: 'include',
      });

      if (!response.ok) {
        if (response.status === 401) {
          console.log('[SessionManager] Session expired during refresh');
          setIsExpired(true);
          return;
        }
        return;
      }

      // Re-check status after refresh
      await checkSessionStatus();
    } catch (error) {
      console.error('[SessionManager] Error refreshing session:', error);
    }
  }, [checkSessionStatus]);

  // Set up periodic status checks
  useEffect(() => {
    // Check if we're on a public route that shouldn't have session management
    const currentPath = window.location.pathname;
    const publicRoutes = ['/login', '/register', '/verify', '/forgot_password', '/reset_password', '/policies'];
    const isPublicRoute = publicRoutes.some(route => currentPath.startsWith(route));
    
    // Only perform session management on protected routes
    if (isPublicRoute) {
      console.log('[SessionManager] On public route, skipping session management');
      return;
    }

    // Load settings and check status immediately
    loadSessionSettings();
    checkSessionStatus();

    // Check status every 30 seconds
    statusCheckIntervalRef.current = window.setInterval(() => {
      checkSessionStatus();
    }, 30000);

    // Update activity periodically based on session settings
    // Use dynamic interval: session timeout minus warning time
    const dynamicActivityInterval = Math.max(
      (sessionSettingsRef.current.timeoutMinutes - sessionSettingsRef.current.warningMinutes) * 60 * 1000,
      300000 // Minimum 5 minutes
    );
    
    activityUpdateIntervalRef.current = window.setInterval(() => {
      console.log('[SessionManager] Periodic activity update (dynamic interval):', {
        interval: Math.round(dynamicActivityInterval / 1000),
        sessionTimeout: sessionSettingsRef.current.timeoutMinutes,
        warningTime: sessionSettingsRef.current.warningMinutes
      });
      updateActivity();
    }, dynamicActivityInterval);

    return () => {
      if (statusCheckIntervalRef.current) {
        clearInterval(statusCheckIntervalRef.current);
      }
      if (activityUpdateIntervalRef.current) {
        clearInterval(activityUpdateIntervalRef.current);
      }
    };
  }, [checkSessionStatus, updateActivity, loadSessionSettings]);

  // Handle session expiration
  useEffect(() => {
    if (isExpired) {
      console.log('[SessionManager] Session expired, logging out');
      setShowWarning(false);
      setTimeRemaining(0);
      // Clear any cached data
      localStorage.clear();
      sessionStorage.clear();
      // Redirect to login page
      window.location.href = '/login';
    }
  }, [isExpired]);

  // Set up activity listeners
  useEffect(() => {
    // Check if we're on a public route that shouldn't have session management
    const currentPath = window.location.pathname;
    const publicRoutes = ['/login', '/register', '/verify', '/forgot_password', '/reset_password', '/policies'];
    const isPublicRoute = publicRoutes.some(route => currentPath.startsWith(route));
    
    // Only set up activity listeners on protected routes
    if (isPublicRoute) {
      console.log('[SessionManager] On public route, skipping activity listeners');
      return;
    }

    const events = ['mousedown', 'scroll', 'touchstart', 'click'];
    
    const handleActivity = () => {
      // Update activity on server when user is active
      updateActivity();
    };

    events.forEach(event => {
      document.addEventListener(event, handleActivity, true);
    });

    return () => {
      events.forEach(event => {
        document.removeEventListener(event, handleActivity, true);
      });
    };
  }, [updateActivity]);

  // Memoize the context value to prevent unnecessary re-renders
  const contextValue = useMemo(() => {
    return {
      showWarning,
      timeRemaining,
      settings: sessionSettings,
      settingsLoaded: true,
      stats: {
        sessionStartTime: 0,
        lastActivityTime: 0,
        lastRefreshTime: 0,
        refreshCount: 0,
        totalSessionDuration: 0,
        timeSinceLastActivity: 0,
        timeSinceLastRefresh: 0,
      },
      updateActivity,
      refreshSession,
      extendSession: refreshSession, // Alias for compatibility
      logoutNow: () => {
        console.log('[SessionManager] Logout requested');
        setIsExpired(true);
        setShowWarning(false);
        setTimeRemaining(0);
        // Clear any cached data
        localStorage.clear();
        sessionStorage.clear();
        // Redirect to login page
        window.location.href = '/login';
      },
    };
  }, [showWarning, timeRemaining, sessionSettings]);

  return (
    <SessionContext.Provider value={contextValue}>
      {children}
    </SessionContext.Provider>
  );
};

export const SessionProvider: React.FC<SessionProviderProps> = ({ 
  children, 
  onSessionExpired 
}) => {
  const onSessionExpiredRef = useRef(onSessionExpired);
  
  useEffect(() => {
    return () => {
      // Cleanup on unmount
    };
  }, []);
  
  // Update the ref when onSessionExpired changes
  useEffect(() => {
    onSessionExpiredRef.current = onSessionExpired;
  }, [onSessionExpired]);

  return (
    <SessionManager>
      {children}
    </SessionManager>
  );
};

export const useSession = (): SessionState => {
  const context = useContext(SessionContext);
  if (context === undefined) {
    throw new Error('useSession must be used within a SessionProvider');
  }
  return context;
};