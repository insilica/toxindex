import React from "react";
import { useNavigate } from "react-router-dom";
import { FaListAlt } from "react-icons/fa";

const Settings: React.FC = () => {
  const navigate = useNavigate();

  return (
    <div className="flex items-center justify-center min-h-screen min-w-full w-full h-screen" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="w-full max-w-2xl p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        <h1 className="text-3xl font-bold mb-8 text-white text-center">Settings</h1>
        <ul className="divide-y divide-gray-700">
          <li className="flex items-center justify-between py-6">
            <span className="text-xl text-white">Environments</span>
            <button
              className="text-white hover:text-purple-200 p-2 rounded-full"
              title="Go to Environments settings"
              onClick={() => navigate("/settings/environments")}
            >
              <FaListAlt size={28} />
            </button>
          </li>
        </ul>
      </div>
    </div>
  );
};

const SettingsGeneral: React.FC = () => {
  return (
    <div className="flex items-center justify-center min-h-screen min-w-full w-full h-screen" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="w-full max-w-2xl p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        <h1 className="text-3xl font-bold mb-8 text-white text-center">General Settings</h1>
        <div className="text-white text-center">This is the General settings page.</div>
      </div>
    </div>
  );
};

export default Settings;
export { SettingsGeneral }; 