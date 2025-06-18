import React, { useState, useEffect } from 'react';
import { useNavigate, Link } from 'react-router-dom';
import { FaUserPlus } from 'react-icons/fa';

const TYPED_TEXT = 'Register';

const RegisterForm: React.FC = () => {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [passwordConfirmation, setPasswordConfirmation] = useState('');
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [typed, setTyped] = useState('');
  const [typing, setTyping] = useState(true);
  const [success, setSuccess] = useState(false);
  const [successMsg, setSuccessMsg] = useState('');
  const [typedSuccess, setTypedSuccess] = useState('');
  const [typingSuccess, setTypingSuccess] = useState(false);
  const navigate = useNavigate();

  useEffect(() => {
    let i = 0;
    setTyped('');
    setTyping(true);
    const interval = setInterval(() => {
      setTyped(TYPED_TEXT.slice(0, i + 1));
      i++;
      if (i === TYPED_TEXT.length) {
        clearInterval(interval);
        setTyping(false);
      }
    }, 70);
    return () => clearInterval(interval);
  }, []);

  useEffect(() => {
    if (success && successMsg) {
      setTypedSuccess('');
      setTypingSuccess(true);
      let i = 0;
      const interval = setInterval(() => {
        setTypedSuccess(successMsg.slice(0, i + 1));
        i++;
        if (i === successMsg.length) {
          clearInterval(interval);
          setTypingSuccess(false);
        }
      }, 40);
      return () => clearInterval(interval);
    }
  }, [success, successMsg]);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);
    setLoading(true);
    if (password !== passwordConfirmation) {
      setError('Passwords do not match.');
      setLoading(false);
      return;
    }
    try {
      const res = await fetch('/api/auth/register', {
        method: 'POST',
        headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
        body: new URLSearchParams({ email, password, password_confirmation: passwordConfirmation }).toString(),
        credentials: 'include',
      });
      if (res.ok) {
        const data = await res.json();
        setSuccess(true);
        setSuccessMsg(data.message || 'Check your email to verify your account.');
        setTimeout(() => navigate('/login'), 5000);
      } else {
        setError('Registration failed.');
      }
    } catch (err) {
      setError('Registration failed.');
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
        .typewriter-cursor {
          display: inline-block;
          width: 1ch;
          color: #6ee7b7;
          font-weight: bold;
          animation: blink 1s steps(1) infinite;
        }
        @keyframes blink {
          0%, 100% { opacity: 1; }
          50% { opacity: 0; }
        }
      `}</style>
      <div
        className="bg-gray-950 rounded-2xl shadow-2xl p-10 w-full max-w-md flex flex-col items-center"
        style={{
          animation: 'fadeInScale 0.7s cubic-bezier(0.4,0,0.2,1) both',
          boxShadow: '0 8px 32px 0 rgba(34,197,94,0.15), 0 1.5px 8px 0 rgba(0,0,0,0.10)'
        }}
      >
        <FaUserPlus className="text-green-400 text-5xl mb-4 drop-shadow-lg animate-bounce-slow" />
        <h2 className="text-3xl font-extrabold mb-2 text-white tracking-tight" style={{ minHeight: 44 }}>
          {typed}
          {typing && <span className="typewriter-cursor">|</span>}
        </h2>
        {success ? (
          <div className="text-green-400 font-bold text-lg mb-2 animate-fade-in" style={{ minHeight: 44 }}>
            {typedSuccess}
            {typingSuccess && <span className="typewriter-cursor">|</span>}
            {!typingSuccess && <div className="text-gray-300 text-sm mt-2 animate-fade-in">
              <br />
              Redirecting to login...</div>}
          </div>
        ) : (
          <form onSubmit={handleSubmit} className="w-full flex flex-col gap-4 mt-2">
            <label className="block text-left text-gray-300 font-semibold">Email</label>
            <input
              type="email"
              className="w-full p-3 rounded-lg border border-gray-700 bg-gray-800 text-white focus:outline-none focus:ring-2 focus:ring-green-400"
              value={email}
              onChange={e => setEmail(e.target.value)}
              required
              disabled={loading}
              placeholder="you@email.com"
            />
            <label className="block text-left text-gray-300 font-semibold">Password</label>
            <input
              type="password"
              className="w-full p-3 rounded-lg border border-gray-700 bg-gray-800 text-white focus:outline-none focus:ring-2 focus:ring-green-400"
              value={password}
              onChange={e => setPassword(e.target.value)}
              required
              disabled={loading}
              placeholder="Create a password"
            />
            <label className="block text-left text-gray-300 font-semibold">Confirm Password</label>
            <input
              type="password"
              className="w-full p-3 rounded-lg border border-gray-700 bg-gray-800 text-white focus:outline-none focus:ring-2 focus:ring-green-400"
              value={passwordConfirmation}
              onChange={e => setPasswordConfirmation(e.target.value)}
              required
              disabled={loading}
              placeholder="Repeat your password"
            />
            {error && <div className="text-red-400 mb-2 text-center animate-fade-in">{error}</div>}
            <button
              type="submit"
              className="w-full bg-gradient-to-r from-green-500 to-green-700 text-white py-3 rounded-lg font-bold shadow-lg hover:from-green-600 hover:to-green-800 transition disabled:opacity-50"
              disabled={loading}
            >
              {loading ? 'Registering...' : 'Register'}
            </button>
          </form>
        )}
        <div className="flex flex-row justify-between w-full mt-6">
          <Link to="/login" className="text-green-400 hover:underline text-sm">Back to Login</Link>
        </div>
      </div>
    </div>
  );
};

export default RegisterForm; 