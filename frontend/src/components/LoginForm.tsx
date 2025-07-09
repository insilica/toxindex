import React, { useState, useEffect } from 'react';
import { useNavigate, Link } from 'react-router-dom';
import { FaLock } from 'react-icons/fa';

const TYPED_TEXT = 'Login';

const LoginForm: React.FC = () => {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [typed, setTyped] = useState('');
  const [typing, setTyping] = useState(true);
  const [errorAnim, setErrorAnim] = useState(false);
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

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);
    setLoading(true);
    try {
      const res = await fetch('/api/auth/login', {
        method: 'POST',
        headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
        body: new URLSearchParams({ email, password }).toString(),
        credentials: 'include',
      });
      if (res.ok) {
        navigate('/');
      } else {
        const data = await res.json();
        setError(data.error || 'Invalid email or password.');
        setErrorAnim(true);
        setTimeout(() => {
          setErrorAnim(false);
          setError(null);
        }, 3000);
      }
    } catch (err) {
      setError('Login failed.');
      setErrorAnim(true);
      setTimeout(() => {
        setErrorAnim(false);
        setError(null);
      }, 3000);
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
        @keyframes shake {
          0% { transform: translateX(0); }
          15% { transform: translateX(-8px); }
          30% { transform: translateX(8px); }
          45% { transform: translateX(-6px); }
          60% { transform: translateX(6px); }
          75% { transform: translateX(-4px); }
          90% { transform: translateX(4px); }
          100% { transform: translateX(0); }
        }
        .shake-glow {
          animation: shake 0.6s cubic-bezier(0.4,0,0.2,1);
          box-shadow: 0 0 24px 6px #fde047, 0 8px 32px 0 rgba(34,197,94,0.15), 0 1.5px 8px 0 rgba(0,0,0,0.10) !important;
        }
        .lock-error {
          color: #f87171 !important;
          filter: drop-shadow(0 0 8px #f87171);
          transition: color 0.2s, filter 0.2s;
        }
        .login-button-normal {
          background: linear-gradient(to right,rgb(242, 234, 214),rgb(183, 193, 207)) !important;
          color:rgb(49, 50, 51) !important;
          border: 2px solid #6b7280 !important;
        }
        .login-error-shake {
          animation: shake 0.6s cubic-bezier(0.4,0,0.2,1) 0s 5;
          background: linear-gradient(90deg, #f87171 0%, #dc2626 100%) !important;
          color: #fff !important;
          border: 2px solid #f87171 !important;
        }
        .enlarge-on-error {
          transform: scale(1.3);
          transition: transform 0.3s cubic-bezier(0.4,0,0.2,1);
        }
      `}</style>
      <div
        className={`bg-gray-950 rounded-2xl shadow-2xl p-10 w-full max-w-md flex flex-col items-center${errorAnim ? ' shake-glow' : ''}`}
        style={{
          animation: 'fadeInScale 0.7s cubic-bezier(0.4,0,0.2,1) both',
          boxShadow: errorAnim
            ? '0 0 24px 6px #fde047, 0 8px 32px 0 rgba(34,197,94,0.15), 0 1.5px 8px 0 rgba(0,0,0,0.10)'
            : '0 8px 32px 0 rgba(34,197,94,0.15), 0 1.5px 8px 0 rgba(0,0,0,0.10)'
        }}
      >
        <FaLock
          className={`text-5xl mb-4 drop-shadow-lg animate-bounce-slow${errorAnim ? ' lock-error' : ' text-sky-200'}`}
        />
        <h2 className="text-3xl font-extrabold mb-2 text-white tracking-tight" style={{ minHeight: 44 }}>
          {typed}
          {typing && <span className="typewriter-cursor">|</span>}
        </h2>
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
            placeholder="Your password"
          />
          {/* Error message placeholder to prevent layout shift */}
          <div style={{ minHeight: 24, marginBottom: 8, textAlign: 'center' }}>
            {error ? (
              <span className="text-red-400 animate-fade-in">{error}</span>
            ) : (
              <span style={{ visibility: 'hidden' }}>placeholder</span>
            )}
          </div>
          <button
            type="submit"
            className={`login-button-normal py-3 px-6 rounded-lg font-bold shadow-lg hover:from-gray-600 hover:to-gray-700 transition-all duration-200 disabled:opacity-50${errorAnim ? ' login-error-shake' : ''}`}
            disabled={loading}
          >
            {loading ? 'Logging in...' : 'Login'}
          </button>
        </form>
        <div className="flex flex-row justify-between w-full mt-6">
          <Link
            to="/register"
            className={`text-neutral-100 hover:underline text-sm font-semibold px-4 py-2 rounded-lg transition hover:bg-green-900/30${errorAnim ? ' enlarge-on-error' : ''}`}
          >
            Create Account
          </Link>
          <Link
            to="/forgot_password"
            className={`text-neutral-100 hover:underline text-sm font-semibold px-4 py-2 rounded-lg transition hover:bg-green-900/30${errorAnim ? ' enlarge-on-error' : ''}`}
          >
            Forgot Password?
          </Link>
        </div>
      </div>
      <footer className="absolute bottom-0 left-0 right-0 p-4 text-center text-neutral-100 text-sm">
        <span className="text-green-400 font-black">TOXINDEX</span> Copyright © Insilica. All Rights Reserved.
      </footer>
    </div>
  );
};

export default LoginForm;