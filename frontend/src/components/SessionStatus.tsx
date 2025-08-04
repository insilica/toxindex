import React, { useState } from 'react';
import { Clock, RefreshCw, AlertTriangle, Activity } from 'lucide-react';
import { useSession } from '../context/SessionContext';

const SessionStatus: React.FC = () => {
  // Get session state from context
  const {
    showWarning,
    timeRemaining,
    timeUntilTimeout,
    settings,
    stats,
    extendSession
  } = useSession();

  const [isExpanded, setIsExpanded] = useState(false);

  const formatDuration = (ms: number): string => {
    const seconds = Math.floor(ms / 1000);
    const minutes = Math.floor(seconds / 60);
    const hours = Math.floor(minutes / 60);
    
    if (hours > 0) {
      return `${hours}h ${minutes % 60}m ${seconds % 60}s`;
    } else if (minutes > 0) {
      return `${minutes}m ${seconds % 60}s`;
    } else {
      return `${seconds}s`;
    }
  };

  const timeoutProgress = (stats.timeSinceLastActivity / (settings.timeoutMinutes * 60 * 1000)) * 100;

  return (
    <div className={`border rounded-lg p-4 ${showWarning ? 'border-orange-500 bg-orange-50' : 'border-gray-200 bg-white'}`}>
      <div className="flex items-center justify-between">
        <div className="flex items-center space-x-2">
          <div className={`p-2 rounded-full ${showWarning ? 'bg-orange-100' : 'bg-blue-100'}`}>
            {showWarning ? (
              <AlertTriangle className={`h-4 w-4 ${showWarning ? 'text-orange-600' : 'text-blue-600'}`} />
            ) : (
              <Clock className="h-4 w-4 text-blue-600" />
            )}
          </div>
          <div>
            <h3 className={`font-medium ${showWarning ? 'text-orange-800' : 'text-gray-900'}`}>
              Session Status
            </h3>
            <p className={`text-sm ${showWarning ? 'text-orange-600' : 'text-gray-600'}`}>
              {showWarning ? (
                `⚠️ Session expires in ${formatDuration(timeRemaining * 1000)}`
              ) : (
                `Active - ${formatDuration(timeUntilTimeout)} remaining`
              )}
            </p>
          </div>
        </div>
        
        <div className="flex items-center space-x-2">
          {showWarning && (
            <button
              onClick={extendSession}
              className="px-3 py-1 bg-orange-600 text-white rounded text-sm hover:bg-orange-700 flex items-center space-x-1"
            >
              <RefreshCw className="h-3 w-3" />
              <span>Extend</span>
            </button>
          )}
          <button
            onClick={() => setIsExpanded(!isExpanded)}
            className="text-gray-400 hover:text-gray-600"
          >
            <Activity className="h-4 w-4" />
          </button>
        </div>
      </div>

      {/* Progress bar */}
      <div className="mt-3">
        <div className="w-full bg-gray-200 rounded-full h-2">
          <div
            className={`h-2 rounded-full transition-all duration-1000 ${
              timeoutProgress > 80 ? 'bg-red-500' : 
              timeoutProgress > 60 ? 'bg-orange-500' : 'bg-green-500'
            }`}
            style={{ width: `${Math.min(100, timeoutProgress)}%` }}
          />
        </div>
      </div>

      {/* Expanded details */}
      {isExpanded && (
        <div className="mt-4 pt-4 border-t border-gray-200 space-y-3">
          <div className="grid grid-cols-2 gap-4 text-sm">
            <div>
              <div className="flex items-center space-x-1 text-gray-600">
                <Clock className="h-3 w-3" />
                <span>Session Duration</span>
              </div>
              <p className="font-medium">{formatDuration(stats.totalSessionDuration)}</p>
            </div>
            
            <div>
              <div className="flex items-center space-x-1 text-gray-600">
                <Activity className="h-3 w-3" />
                <span>Last Activity</span>
              </div>
              <p className="font-medium">{formatDuration(stats.timeSinceLastActivity)} ago</p>
            </div>
            
            <div>
              <div className="flex items-center space-x-1 text-gray-600">
                <RefreshCw className="h-3 w-3" />
                <span>Last Refresh</span>
              </div>
              <p className="font-medium">{formatDuration(stats.timeSinceLastRefresh)} ago</p>
            </div>
            
            <div>
              <div className="flex items-center space-x-1 text-gray-600">
                <RefreshCw className="h-3 w-3" />
                <span>Refresh Count</span>
              </div>
              <p className="font-medium">{stats.refreshCount} times</p>
            </div>
          </div>
          
          <div className="pt-2 border-t border-gray-100">
            <div className="text-xs text-gray-500 space-y-1">
              <p>Timeout: {settings.timeoutMinutes} minutes</p>
              <p>Warning: {settings.warningMinutes} minutes before timeout</p>
              <p>Auto-refresh: Every {settings.refreshIntervalMinutes} minutes</p>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default SessionStatus;