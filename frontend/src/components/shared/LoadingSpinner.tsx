import React from 'react';
import { FaSpinner } from 'react-icons/fa';

interface LoadingSpinnerProps {
  size?: 'small' | 'medium' | 'large';
  className?: string;
  text?: string;
  showTimer?: boolean;
  startTime?: number; // UNIX timestamp in seconds or ms
  workflowId?: number;
  status?: string;
}

// Status-to-step mapping for each workflow/tool
const STATUS_STEP_MAP: Record<number, Record<string, number>> = {
  1: { // probra
    "starting": 1,
    "checking cache": 2,
    "using cache": 2,
    "running agent": 2,
    "agent complete": 3,
    "sending message": 3,
    "file created": 4,
    "done": 4,
    "error": 0,
  },
  2: { // plain_openai_task
    "checking cache": 1,
    "using cache": 2,
    "calling openai": 2,
    "openai complete": 3,
    "publishing message": 4,
    "done": 4,
    "error": 0,
  },
  3: { // openai_json_schema_task
    "checking cache": 1,
    "using cache": 2,
    "calling openai": 2,
    "openai complete": 3,
    "publishing message": 4,
    "done": 4,
    "error": 0,
  },
  4: { // raptool_task
    "fetching file": 1,
    "converting csv": 2,
    "parsing chemicals": 3,
    "publishing result": 4,
    "done": 4,
    "error": 0,
  },
};

const ProgressIndicator = ({ step, total = 4 }: { step: number, total?: number }) => (
  <span className="ml-2 flex items-center">
    {[...Array(total)].map((_, i) => (
      <span
        key={i}
        className={`mx-0.5 w-2 h-2 rounded-full ${i < step ? 'bg-green-400' : 'bg-gray-700'} inline-block`}
        style={{ border: '1px solid #6ee7b7' }}
      />
    ))}
  </span>
);

const LoadingSpinner: React.FC<LoadingSpinnerProps> = ({
  size = 'medium',
  className = '',
  text = '',
  showTimer = false,
  startTime,
  workflowId,
  status
}) => {
  const sizeClass = {
    small: 'text-xl',
    medium: 'text-3xl',
    large: 'text-4xl'
  }[size];

  // If startTime is provided, calculate elapsed seconds from it
  const getInitialSeconds = () => {
    if (!startTime) return 0;
    const now = Date.now();
    // If startTime is in seconds, convert to ms
    const startMs = startTime > 1e12 ? startTime : startTime * 1000;
    return Math.floor((now - startMs) / 1000);
  };

  const [seconds, setSeconds] = React.useState(getInitialSeconds());
  React.useEffect(() => {
    if (!showTimer) return;
    const update = () => setSeconds(getInitialSeconds());
    update(); // set immediately
    const interval = setInterval(update, 1000);
    return () => clearInterval(interval);
  }, [showTimer, startTime]);

  // Determine progress step
  const step = workflowId && status
    ? STATUS_STEP_MAP[workflowId]?.[status] ?? 0
    : 0;

  return (
    <div className={`flex items-center justify-center flex-col gap-3 ${className}`}>
      <div className="flex items-center gap-2">
        <FaSpinner className={`animate-spin ${sizeClass} text-gray-500`} />
        {showTimer && (
          <span className="text-xs text-gray-400 font-mono" aria-label={`It's been ${seconds} seconds`}>
            {seconds}secs
          </span>
        )}
        {/* Progress indicator to the right of spinner/timer */}
        {workflowId && status && <ProgressIndicator step={step} />}
      </div>
      {text && <span className="text-gray-400">{text}</span>}
    </div>
  );
};

export default LoadingSpinner; 