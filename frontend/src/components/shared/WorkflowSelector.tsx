import React, { useEffect, useRef, useState } from 'react';
import { useModel } from '../../context/ModelContext';
import { FaHammer } from 'react-icons/fa';
import { WORKFLOWS, getWorkflowByFrontendId, loadWorkflows } from './workflows';

interface WorkflowSelectorProps {
  className?: string;
}

const WorkflowSelector: React.FC<WorkflowSelectorProps> = ({ className = '' }) => {
  const selectRef = useRef<HTMLSelectElement>(null);
  const spanRef = useRef<HTMLSpanElement>(null);
  const [selectWidth, setSelectWidth] = useState(120);
  const [focused, setFocused] = useState(false);
  const [workflows, setWorkflows] = useState(WORKFLOWS);
  const { selectedModel, setSelectedModel } = useModel();

  // Handle clicking outside to clear focus
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (selectRef.current && !selectRef.current.contains(event.target as Node)) {
        // Force blur and clear focus state
        selectRef.current.blur();
        setFocused(false);
        // Additional timeout to ensure state is cleared
        setTimeout(() => {
          setFocused(false);
        }, 10);
      }
    };

    document.addEventListener('click', handleClickOutside);
    return () => {
      document.removeEventListener('click', handleClickOutside);
    };
  }, []);

  useEffect(() => {
    // Force refresh workflows on component mount
    loadWorkflows(true).then(loadedWorkflows => {
      setWorkflows(loadedWorkflows);
    });
  }, []);

  useEffect(() => {
    if (spanRef.current) {
      const measuredWidth = spanRef.current.offsetWidth;
      // Ensure minimum width of 120px, and add padding
      const calculatedWidth = Math.max(measuredWidth + 40, 200);
      setSelectWidth(calculatedWidth);
    }
  }, [selectedModel, workflows]);

  // Ensure selectedModel is valid for available workflows
  useEffect(() => {
    if (workflows.length > 0) {
      const availableWorkflowIds = workflows.map(w => w.frontend_id);
      if (!availableWorkflowIds.includes(selectedModel)) {
        // If current selection is not available, select the first available workflow
        setSelectedModel(availableWorkflowIds[0]);
      }
    }
  }, [workflows, selectedModel, setSelectedModel]);

  if (workflows.length === 0) {
    return (
      <div className="relative group flex items-center gap-0" style={{ minWidth: 200, maxWidth: 300 }}>
        <FaHammer className="text-blue-400" style={{ fontSize: '1.6rem', flexShrink: 0 }} title="Workflow" />
        <div className="text-white text-sm">Loading...</div>
      </div>
    );
  }

  return (
    <div className="relative group flex items-center gap-0" style={{ minWidth: 200, maxWidth: 300 }}>
      <FaHammer className="text-blue-400" style={{ fontSize: '1.6rem', flexShrink: 0 }} title="Workflow" />
      {/* Hidden span to measure width */}
      <span
        ref={spanRef}
        style={{
          position: 'absolute',
          visibility: 'hidden',
          whiteSpace: 'nowrap',
          fontWeight: 'bold',
          fontSize: '0.875rem',
          fontFamily: 'inherit',
          padding: '0 16px',
        }}
      >
        {getWorkflowByFrontendId(selectedModel)?.label || workflows[0]?.label || 'Select Workflow'}
      </span>
      <select
        ref={selectRef}
        className={`font-bold text-white text-base leading-tight pl-2 pr-4 py-2 rounded-full appearance-none bg-transparent border-none focus:outline-none transition-all duration-150 group-hover:bg-black group-hover:bg-opacity-60 group-hover:border group-hover:border-gray-700 group-hover:pl-4 group-hover:pr-12 group-hover:cursor-pointer ${focused ? 'bg-black bg-opacity-60 border border-gray-700 pl-4 pr-12' : ''} w-full ${className}`}
        style={{
          width: `${selectWidth}px`,
          minWidth: '50px',
          maxWidth: '260px',
          height: '2.5rem',
          borderRadius: '999px',
          position: 'relative',
          zIndex: 1,
          backgroundColor: 'transparent',
          cursor: 'pointer',
          paddingRight: '2.2rem',
          lineHeight: '1.5rem',
        }}
        value={selectedModel}
        onChange={e => setSelectedModel(e.target.value)}
        onFocus={() => setFocused(true)}
      >
        {workflows.map(w => (
          <option key={w.frontend_id} value={w.frontend_id} className="!bg-gray-800 !text-white" style={{ paddingLeft: '1rem' }}>
            {w.label}
          </option>
        ))}
      </select>
      <span
        className="pointer-events-none absolute right-2 top-1/2 transform -translate-y-1/2 text-gray-400 group-hover:opacity-100 group-focus-within:opacity-100 opacity-0 transition-opacity duration-150"
        style={{ right: '1.2rem' }}
      >
        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round">
          <path d="M7 10l5 5 5-5" />
        </svg>
      </span>
    </div>
  );
};

export default WorkflowSelector; 