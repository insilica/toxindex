import React, { useState } from 'react';
import { useParams, useNavigate, Link } from 'react-router-dom';
import { FaLock } from 'react-icons/fa';

const ResetPasswordPage: React.FC = () => {
  const { token } = useParams<{ token: string }>();
  const [password, setPassword] = useState('');
  const [success, setSuccess] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const navigate = useNavigate();

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);
    setLoading(true);
    try {
      const res = await fetch(`/api/auth/reset_password/${token}`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ password }),
      });
      if (res.ok) {
        const data = await res.json();
        if (data.success) {
          setSuccess(true);
          setTimeout(() => navigate('/login'), 2000);
        } else {
          setError(data.error || 'Failed to reset password.');
        }
      } else {
        let errorMsg = 'Failed to reset password.';
        try {
          const data = await res.json();
          errorMsg = data.error || errorMsg;
        } catch {
          const text = await res.text();
          if (text) errorMsg = text;
        }
        setError(errorMsg);
      }
    } catch (err) {
      setError('Failed to reset password.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen flex items-center justify-center bg-gradient-to-br from-green-900 via-gray-900 to-black">
      <style>{`
        @keyframes fadeInScale {
          0% { opacity: 0; transform: scale(0.95); }
          100% { opacity: 1; transform: scale(1); }
        }
      `}</style>
      <div
        className="bg-gray-950 rounded-2xl shadow-2xl p-10 w-full max-w-md flex flex-col items-center"
        style={{
          animation: 'fadeInScale 0.7s cubic-bezier(0.4,0,0.2,1) both',
          boxShadow: '0 8px 32px 0 rgba(34,197,94,0.15), 0 1.5px 8px 0 rgba(0,0,0,0.10)'
        }}
      >
        <FaLock className="text-green-400 text-5xl mb-4 drop-shadow-lg animate-bounce-slow" />
        <h2 className="text-3xl font-extrabold mb-2 text-white tracking-tight">Reset Password</h2>
        <p className="text-gray-400 mb-6 text-center">Enter your new password below.</p>
        {success ? (
          <div className="text-green-400 mb-4 text-center animate-fade-in">Password reset! Redirecting to login...</div>
        ) : (
          <form onSubmit={handleSubmit} className="w-full flex flex-col gap-4">
            <label className="block text-left text-gray-300 font-semibold">New Password</label>
            <input
              type="password"
              className="w-full p-3 rounded-lg border border-gray-700 bg-gray-800 text-white focus:outline-none focus:ring-2 focus:ring-green-400"
              value={password}
              onChange={e => setPassword(e.target.value)}
              required
              disabled={loading}
              placeholder="Enter new password"
            />
            {error && <div className="text-red-400 mb-2 text-center animate-fade-in">{error}</div>}
            <button
              type="submit"
              className="w-full !bg-gray-800 text-white py-3 rounded-lg font-bold shadow-lg transition disabled:opacity-50"
              disabled={loading}
            >
              {loading ? 'Resetting...' : 'Reset Password'}
            </button>
          </form>
        )}
        <Link to="/login" className="mt-6 text-green-400 hover:underline text-sm">Back to Login</Link>
      </div>
      <style>{`
        .animate-bounce-slow {
          animation: bounce 2.2s infinite cubic-bezier(0.4,0,0.2,1);
        }
        @keyframes bounce {
          0%, 100% { transform: translateY(0); }
          50% { transform: translateY(-10px); }
        }
        .animate-fade-in {
          animation: fadeInScale 0.7s cubic-bezier(0.4,0,0.2,1) both;
        }
      `}</style>
    </div>
  );
};

export default ResetPasswordPage; 