import React from "react";
import { useNavigate } from "react-router-dom";

const VerifyPage: React.FC = () => {
  const navigate = useNavigate();
  return (
    <div className="min-h-screen flex items-center justify-center" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)' }}>
      <div className="max-w-md w-full p-8 bg-white rounded-2xl shadow-2xl flex flex-col items-center" style={{ boxShadow: '0 8px 32px 0 rgba(34,197,94,0.10)' }}>
        <div className="flex items-center mb-4">
          <svg width="36" height="36" viewBox="0 0 24 24" fill="none" stroke="#22c55e" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
            <circle cx="12" cy="12" r="10" />
            <path d="M9 12l2 2l4-4" />
          </svg>
          <h2 className="text-2xl font-bold text-gray-900">Verify Your Email</h2>
        </div>
        <p className="text-gray-700 text-center mb-6">
          Please check your email and click the verification link to activate your account.
        </p>
        <button
          onClick={() => navigate('/login')}
          className="mt-2 px-6 py-2 bg-green-600 hover:bg-green-700 text-white rounded-full font-semibold shadow transition"
        >
          Go to Login
        </button>
      </div>
    </div>
  );
};

export default VerifyPage; 