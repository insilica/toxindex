import React, { useEffect, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';

const VerifyPage: React.FC = () => {
  console.log("VerifyPage is being rendered");
  const { token } = useParams<{ token: string }>();
  const [status, setStatus] = useState<'pending' | 'success' | 'error'>('pending');
  const [message, setMessage] = useState('');
  const navigate = useNavigate();

  useEffect(() => {
    if (!token) {
      setStatus('error');
      setMessage('No verification token provided.');
      return;
    }

    fetch(`/api/auth/verification/${token}`)
      .then(async res => {
        const data = await res.json();
        if (res.ok) {
          setStatus('success');
          setMessage(data.message || 'Your account has been activated!');
          setTimeout(() => navigate('/login'), 3000);
        } else {
          setStatus('error');
          setMessage(data.error || 'Verification failed.');
        }
      })
      .catch(error => {
        console.error('Verification error:', error);
        setStatus('error');
        setMessage('An error occurred during verification. Please try again.');
      });
  }, [token, navigate]);

  return (
    <div className="min-h-screen flex items-center justify-center bg-gradient-to-br from-green-900 via-gray-900 to-black">
      <div className="bg-gray-950 rounded-2xl shadow-2xl p-10 w-full max-w-md flex flex-col items-center">
        <h2 className="text-3xl font-extrabold mb-4 text-white tracking-tight">Account Verification</h2>
        {status === 'pending' && <div className="text-gray-300">Verifying your account...</div>}
        {status === 'success' && <div className="text-green-400 font-bold text-lg mb-2">{message} Redirecting to login...</div>}
        {status === 'error' && (
          <div className="text-red-400 font-bold text-lg mb-2">
            {message}
            <div className="mt-4">
              <button 
                onClick={() => navigate('/login')} 
                className="px-4 py-2 bg-gray-800 hover:bg-gray-700 text-white rounded-lg transition-colors"
              >
                Return to Login
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default VerifyPage; 