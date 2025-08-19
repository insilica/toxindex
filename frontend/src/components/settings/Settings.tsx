import React, { useState, useRef } from "react";
import { FaListAlt } from "react-icons/fa";
import HomeButton from '../shared/HomeButton';

const DEFAULT_INSTRUCTIONS = "Example: Run tests and linters for every code change but not when changing code comments or documentation";

const Settings: React.FC = () => {
  const [instructions, setInstructions] = useState<string>("");
  const [focused, setFocused] = useState(false);
  const textareaRef = useRef<HTMLTextAreaElement>(null);

  return (
    <div className="max-w-6xl mx-auto p-6 relative" style={{ paddingLeft: '8rem' }}>
      <div className="w-full max-w-2xl p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        <h1 className="text-xs font-medium text-gray-200 text-center mb-2 tracking-wide uppercase">Settings</h1>
        <div className="mb-8">
          <h2 className="text-xl font-semibold text-white mb-2 flex items-center gap-2">
            <FaListAlt className="text-green-400" /> Custom Instructions
          </h2>
          <p className="text-gray-300 mb-4 text-sm">
            Set custom instructions for how automated tasks should behave. This can help tailor the workflow to your team's needs.
          </p>
          <div className="relative w-full">
            <textarea
              ref={textareaRef}
              className="w-full min-h-[80px] max-h-[200px] p-3 rounded bg-gray-900 bg-opacity-60 text-white border border-gray-700 focus:border-green-400 focus:outline-none text-base resize-vertical shadow"
              value={instructions}
              onChange={e => setInstructions(e.target.value)}
              onFocus={() => setFocused(true)}
              onBlur={() => setFocused(false)}
            />
            {(!instructions && !focused) && (
              <div
                className="absolute left-3 top-3 text-gray-500 pointer-events-none select-none"
                style={{ fontSize: "1rem" }}
              >
                {DEFAULT_INSTRUCTIONS}
              </div>
            )}
          </div>
        </div>
        <div className="text-gray-400 text-sm text-center mt-8">
          Use the sidebar to access more settings, such as Environments and Data Controls.
        </div>
      </div>
      
      {/* Home button positioned relative to settings content */}
      <div 
        className="absolute transition-all duration-300 z-50"
        style={{
          left: '2rem',
          top: '1.5rem',
          border: 'none',
          padding: 0
        }}
      >
        <HomeButton
          color="#16a34a"
          hoverColor="#2563eb"
          aria-label="Go back"
        />
      </div>
    </div>
  );
};

export default Settings;
