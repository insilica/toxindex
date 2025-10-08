import React, { useEffect, useState, useRef } from 'react';
import { useParams } from 'react-router-dom';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { useEnvironment } from "../context/EnvironmentContext";
import { useModel } from "../context/ModelContext";
import ChatInputBar from './shared/ChatInputBar';
import HomeButton from './shared/HomeButton';
import LoadingSpinner from './shared/LoadingSpinner';
import { getWorkflowId } from './shared/workflows';
import { useSocket } from '../context/SocketContext';

let mountCount = 0;

const ChatSession = () => {
  const { sessionId } = useParams<{ sessionId: string }>();
  const [messages, setMessages] = useState<any[]>([]);
  const [loading, setLoading] = useState(true);
  const [sending, setSending] = useState(false);
  const [input, setInput] = useState('');
  const [error, setError] = useState<string | null>(null);
  const messagesEndRef = useRef<HTMLDivElement>(null);
  const { selectedEnv } = useEnvironment();
  const { selectedModel } = useModel();
  const [fileId, setFileId] = useState<string | undefined>(undefined);
  const [fileName, setFileName] = useState<string | undefined>(undefined);
  const { socket, isConnected, connect } = useSocket();
  const [chatTitle, setChatTitle] = useState<string>('Chat Session');
  const [chatCreatedAt, setChatCreatedAt] = useState<string | null>(null);
  const getMessageKey = (m: any): string => {
    if (!m) return 'null';
    return m.message_id || `${m.role || 'unknown'}::${m.content || ''}`;
  };

  // Keep only first assistant message after each user message
  const dedupeMessages = (messages: any[]): any[] => {
    const result: any[] = [];
    let lastUserIndex = -1;
    
    for (let i = 0; i < messages.length; i++) {
      const msg = messages[i];
      result.push(msg);
      
      if (msg.role === 'user') {
        lastUserIndex = i;
      } else if (msg.role === 'assistant' && lastUserIndex >= 0) {
        // This is an assistant message after a user message
        // Remove any subsequent assistant messages until next user message
        let j = i + 1;
        while (j < messages.length && messages[j].role === 'assistant') {
          j++; // Skip duplicate assistant messages
        }
        i = j - 1; // Adjust loop index
      }
    }
    
    return result;
  };
  const isDev = (import.meta as any)?.env?.MODE !== 'production';

  console.log("ChatSession mounted", sessionId);

  // Connect Socket.IO once on mount
  useEffect(() => {
    console.log('[ChatSession] Connecting Socket.IO...');
    connect();
  }, []); // Empty dependency array - only run once

  // Debug socket connection state
  useEffect(() => {
    console.log('[ChatSession] Socket state changed:', { 
      hasSocket: !!socket, 
      isConnected, 
      sessionId 
    });
  }, [socket, isConnected, sessionId]);

  useEffect(() => {
    mountCount += 1;
    console.log("ChatSession useEffect mount count:", mountCount);
  }, []);

  useEffect(() => {
    if (!socket || !sessionId) return;

    // Join room if already connected
    if (isConnected) {
      console.log('[SocketIO] Already connected, joining room immediately');
      socket.emit('join_chat_session', { session_id: sessionId });
    }

    // Re-join room on every connect (handles reconnects)
    const onConnect = () => {
      console.log('[SocketIO] Connected, joining chat session room:', sessionId);
      socket.emit('join_chat_session', { session_id: sessionId });
    };

    // Handler for new_message
    const handler = (msg: any) => {
      console.log('[SocketIO] Received new_message:', msg);
      console.log('[SocketIO] Message details:', { role: msg?.role, content: msg?.content?.substring(0, 50) + '...' });
      // Append and dedupe to keep only first assistant message after each user message
      setMessages(prev => {
        const newMessages = [...prev, msg];
        return dedupeMessages(newMessages);
      });
    };

    // Listen for room join confirmation
    const onJoinedRoom = (data: any) => {
      console.log('[SocketIO] Successfully joined chat session:', data);
    };

    // Test: listen for ANY socket event to see if socket is working
    const onAnyEvent = (...args: any[]) => {
      console.log('[SocketIO] Received ANY event:', args);
    };
    socket.onAny(onAnyEvent);

    socket.on('connect', onConnect);
    socket.on('new_message', handler);
    socket.on('joined_chat_session', onJoinedRoom);

    // Cleanup handlers
    return () => {
      socket.offAny(onAnyEvent);
      socket.off('connect', onConnect);
      socket.off('new_message', handler);
      socket.off('joined_chat_session', onJoinedRoom);
    };
  }, [sessionId, socket, isConnected]);

  useEffect(() => {
    if (!sessionId) return;
    setLoading(true);
    fetch(`/api/chat_sessions/${sessionId}/messages`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        const list = Array.isArray(data.messages) ? data.messages : [];
        // Keep only first assistant message after each user message
        const deduped = dedupeMessages(list);
        setMessages(deduped);
      })
      .catch(err => {
        setMessages([]);
        setError("Failed to load messages.");
        console.error("Error loading messages:", err);
      })
      .finally(() => setLoading(false));
  }, [sessionId]);

  // Fetch chat session metadata (title and created_at)
  useEffect(() => {
    if (!sessionId) return;
    console.log('[ChatSession] Fetching session metadata for:', sessionId);
    fetch(`/api/chat_sessions/${sessionId}`, { credentials: 'include' })
      .then(res => {
        console.log('[ChatSession] Session metadata response:', res.status, res.ok);
        if (!res.ok) {
          console.error('[ChatSession] API call failed:', res.status, res.statusText);
          return null;
        }
        return res.json();
      })
      .then(data => {
        console.log('[ChatSession] Session metadata data:', data);
        if (!data) {
          console.log('[ChatSession] No data received, keeping defaults');
          return;
        }
        if (data.title) {
          console.log('[ChatSession] Setting title:', data.title);
          setChatTitle(data.title);
        } else {
          console.log('[ChatSession] No title in data, keeping default');
        }
        if (data.created_at) {
          console.log('[ChatSession] Setting created_at:', data.created_at);
          setChatCreatedAt(data.created_at);
        } else {
          console.log('[ChatSession] No created_at in data, keeping null');
        }
      })
      .catch((err) => {
        console.error('[ChatSession] Failed to fetch session metadata:', err);
        // Silently ignore; fallback values already set
      });
  }, [sessionId]);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
  };

  useEffect(() => {
    scrollToBottom();
  }, [messages]);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!input.trim() || sending) return;
    const workflow_id = getWorkflowId(selectedModel);
    // Block submit if workflow 4 and no file selected
    if (workflow_id === 4 && !fileId) {
      setError("You must select a file for this workflow.");
      return;
    }
    setSending(true);
    setError(null);
    if (workflow_id === 0) {
      setError("This workflow is not yet supported.");
      setSending(false);
      return;
    }
    try {
      const res = await fetch('/api/tasks', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          message: input.trim(),
          workflow: workflow_id,
          environment_id: selectedEnv || undefined,
          sid: sessionId,
          file_id: fileId,
        }),
        credentials: 'include'
      });
      const data = await res.json();
      if (data.error) {
        setError(data.error);
      } else {
        setMessages(prev => [...prev, { content: input.trim(), role: 'user' }]);
        setInput('');
        setFileId(undefined);
        setFileName(undefined);
        
        // Fallback: poll for new messages after 3 seconds if no socket message arrives
        // TODO: Remove this fallback once socket events are working reliably
        // setTimeout(() => {
        //   console.log('[ChatSession] Polling for new messages as fallback...');
        //   fetch(`/api/chat_sessions/${sessionId}/messages`, { credentials: 'include' })
        //     .then(res => res.json())
        //     .then(data => {
        //       const list = Array.isArray(data.messages) ? data.messages : [];
        //       const deduped = dedupeMessages(list);
        //       console.log('[ChatSession] Fallback poll completed, messages:', deduped.length);
        //       console.log('[ChatSession] Last message:', deduped[deduped.length - 1]);
        //       setMessages(deduped);
        //     })
        //     .catch(err => console.error('[ChatSession] Fallback poll failed:', err));
        // }, 3000);
      }
    } catch (error) {
      setError('Error sending message.');
      console.error('Error sending message:', error);
    } finally {
      setSending(false);
    }
  };

  return (
    <div className="min-h-screen flex flex-col flex-1 relative" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)' }}>
      {/* Header */}
      <div className="flex items-center gap-4 px-4 py-4 border-b-2 border-gray-600" style={{ marginLeft: 'auto', marginRight: 'auto', width: '100%', maxWidth: '64rem' }}>
        <HomeButton />
        <span className="text-neutral-600 text-xl font-light mx-2">|</span>
        <span className="text-lg font-semibold text-white truncate">{chatTitle}</span>
        <span className="text-neutral-400 ml-4 text-sm whitespace-nowrap">
          {chatCreatedAt ? new Date(chatCreatedAt).toLocaleString() : 'Unknown date'}
        </span>
      </div>
      <div className="flex-1 flex flex-col max-w-5xl mx-auto w-full px-4 pt-6 py-0">
        {/* Chat history box */}
        <div
          className="flex-1 overflow-y-auto w-full max-w-5xl mx-auto px-20 py-0 bg-transparent rounded-lg shadow-inner"
          style={{
            scrollBehavior: 'smooth',
            background: 'transparent',
            minHeight: 0, // allow flex child to shrink
            maxHeight: '75vh', // limit height so only this box scrolls
          }}
        >
          {loading ? (
            <LoadingSpinner text="Loading messages..." />
          ) : (
            <>
              {messages.length === 0 && (
                <div className="text-gray-400 text-center mt-12">No messages yet. Start the conversation!</div>
              )}
              {messages.map((msg, idx) => (
                <div key={idx} className={`flex mb-4 ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`} style={{ width: '100%' }}>
                  {msg.role === 'user' ? (
                    <div
                      className="rounded-2xl px-5 py-3 shadow-md max-w-[70%] bg-gradient-to-br from-purple-600 to-green-500 text-white ml-auto"
                      style={{
                        textAlign: 'right',
                        fontSize: '1rem',
                        wordBreak: 'break-word',
                        borderBottomRightRadius: 8,
                        borderBottomLeftRadius: 24,
                        borderTopLeftRadius: 24,
                        borderTopRightRadius: 24,
                        marginLeft: 32,
                        marginRight: 0,
                        boxShadow: '0 2px 8px 0 rgba(34,197,94,0.10)'
                      }}
                    >
                      <ReactMarkdown remarkPlugins={[remarkGfm]}>{msg.content}</ReactMarkdown>
                      {isDev && (
                        <div className="text-xs text-gray-200 mt-1 opacity-75 break-all" title="dedupe key">
                          key: {getMessageKey(msg)}
                        </div>
                      )}
                    </div>
                  ) : (
                    <div
                      className="text-gray-100 w-full"
                      style={{
                        textAlign: 'left',
                        fontSize: '1rem',
                        wordBreak: 'break-word',
                        marginLeft: 0,
                        lineHeight: 1.8,
                      }}
                    >
                      <ReactMarkdown
                        remarkPlugins={[remarkGfm]}
                        components={{
                          td: ({node, ...props}) => <td className="markdown-table-cell" {...props} />
                        }}
                      >
                        {msg.content}
                      </ReactMarkdown>
                      {isDev && (
                        <div className="text-xs text-gray-500 mt-1 opacity-75 break-all" title="dedupe key">
                          key: {getMessageKey(msg)}
                        </div>
                      )}
                    </div>
                  )}
                </div>
              ))}
              {messages.length > 0 && messages[messages.length-1]?.role === 'user' && (
                <div className="flex justify-start mb-2">
                  <div className="rounded-xl px-4 py-2 bg-gray-700 text-white max-w-[70%] text-sm opacity-80">
                    Task is running...
                  </div>
                </div>
              )}
              <div ref={messagesEndRef} />
            </>
          )}
        </div>
        {/* Selected file display above input bar, left-aligned */}
        <div className="w-full flex flex-col items-center pb-8 relative">
          <div style={{ width: '100%', maxWidth: 800 }}>
            {fileId && (
              <div className="bg-gray-900/10 text-green-200 px-3 py-0 rounded-full text-sm font-medium flex items-center gap-1 mb-2" style={{ minHeight: 40 }}>
                <span className="font-mono font-semibold text-base z-13">Selected file:</span>
                <span className="font-mono truncate max-w-xs text-base z-13">{fileName ? (fileName.length > 24 ? fileName.slice(0, 20) + '...' : fileName) : fileId.slice(0, 10) + '...'}</span>
                <button onClick={() => { setFileId(undefined); setFileName(undefined); }} className="ml-0 text-red-300 hover:text-red-500 z-13" title="Clear file">×</button>
              </div>
            )}
            <ChatInputBar
              value={input}
              onChange={setInput}
              onSubmit={handleSubmit}
              uploading={sending}
              error={error}
              onFilePick={(id: string, name?: string) => {
                setFileId(id);
                setFileName(name);
              }}
            />
          </div>
        </div>
      </div>
    </div>
  );
};

export default ChatSession; 